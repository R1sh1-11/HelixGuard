import pandas as pd
import pytest
from unittest.mock import patch, MagicMock
from src.sanitize import validate_and_clean, sanitize, MISSING_GENOTYPES, VALID_GENOTYPES


# ------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------

def make_df(rows):
    """Helper: build a minimal genome DataFrame like parse_genome() returns."""
    return pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "genotype"])


FAKE_BLOCKLIST = {
    "rs123": "AA",  # flagged SNP
    "rs456": "GG",  # flagged SNP that is also an LD neighbor of rs123
}

FAKE_LD_NEIGHBORS = {
    "rs123": ["rs456", "rs789"],  # rs456 is both blocklisted and a neighbor
    "rs456": ["rs999"],
}


# ------------------------------------------------------------------
# Edge case 1: SNP is both blocklisted AND an LD neighbor
# ------------------------------------------------------------------

class TestBlocklistLDOverlap:
    def test_blocklist_takes_priority_over_ld(self):
        """rs456 is blocklisted AND an LD neighbor of rs123.
        It should get blocklist replacement, NOT be tagged as ld_neighbor."""
        df = make_df([
            {"rsid": "rs123", "chrom": "1", "pos": 100, "genotype": "CT"},
            {"rsid": "rs456", "chrom": "1", "pos": 200, "genotype": "AC"},
            {"rsid": "rs789", "chrom": "1", "pos": 300, "genotype": "TT"},
        ])

        with patch("src.sanitize.load_blocklist", return_value=FAKE_BLOCKLIST), \
             patch("src.sanitize.parse_genome", return_value=df), \
             patch("src.sanitize.get_ld_neighbors", side_effect=lambda rsid: FAKE_LD_NEIGHBORS.get(rsid, [])):

            result = sanitize("fake_path.txt")

        # rs456 should have blocklist genotype "GG", not be ld_neighbor
        rs456_row = result[result["rsid"] == "rs456"].iloc[0]
        assert rs456_row["genotype"] == "GG", "blocklist replacement should apply to rs456"
        assert rs456_row["ld_neighbor"] == False, "rs456 should NOT be tagged as ld_neighbor"

    def test_ld_only_neighbor_is_tagged(self):
        """rs789 is ONLY an LD neighbor, not blocklisted -- should be tagged."""
        df = make_df([
            {"rsid": "rs123", "chrom": "1", "pos": 100, "genotype": "CT"},
            {"rsid": "rs456", "chrom": "1", "pos": 200, "genotype": "AC"},
            {"rsid": "rs789", "chrom": "1", "pos": 300, "genotype": "TT"},
        ])

        with patch("src.sanitize.load_blocklist", return_value=FAKE_BLOCKLIST), \
             patch("src.sanitize.parse_genome", return_value=df), \
             patch("src.sanitize.get_ld_neighbors", side_effect=lambda rsid: FAKE_LD_NEIGHBORS.get(rsid, [])):

            result = sanitize("fake_path.txt")

        rs789_row = result[result["rsid"] == "rs789"].iloc[0]
        assert rs789_row["ld_neighbor"] == True, "rs789 should be tagged as ld_neighbor"
        assert rs789_row["genotype"] == "TT", "rs789 genotype should be unchanged"


# ------------------------------------------------------------------
# Edge case 2: rsIDs not in blocklist pass through unchanged
# ------------------------------------------------------------------

class TestNonBlocklistedPassthrough:
    def test_unflagged_snps_unchanged(self):
        """SNPs not in the blocklist should have identical genotype in output."""
        df = make_df([
            {"rsid": "rs_unflagged_1", "chrom": "2", "pos": 500, "genotype": "AG"},
            {"rsid": "rs_unflagged_2", "chrom": "2", "pos": 600, "genotype": "CC"},
            {"rsid": "rs_unflagged_3", "chrom": "2", "pos": 700, "genotype": "TT"},
        ])

        with patch("src.sanitize.load_blocklist", return_value=FAKE_BLOCKLIST), \
             patch("src.sanitize.parse_genome", return_value=df), \
             patch("src.sanitize.get_ld_neighbors", return_value=[]):

            result = sanitize("fake_path.txt")

        assert result[result["rsid"] == "rs_unflagged_1"].iloc[0]["genotype"] == "AG"
        assert result[result["rsid"] == "rs_unflagged_2"].iloc[0]["genotype"] == "CC"
        assert result[result["rsid"] == "rs_unflagged_3"].iloc[0]["genotype"] == "TT"

    def test_unflagged_snps_not_ld_tagged(self):
        """Unflagged SNPs with no LD relationship should have ld_neighbor=False."""
        df = make_df([
            {"rsid": "rs_unflagged_1", "chrom": "2", "pos": 500, "genotype": "AG"},
        ])

        with patch("src.sanitize.load_blocklist", return_value=FAKE_BLOCKLIST), \
             patch("src.sanitize.parse_genome", return_value=df), \
             patch("src.sanitize.get_ld_neighbors", return_value=[]):

            result = sanitize("fake_path.txt")

        assert result.iloc[0]["ld_neighbor"] == False


# ------------------------------------------------------------------
# Edge case 3: Missing/malformed genotypes don't crash the pipeline
# ------------------------------------------------------------------

class TestMalformedGenotypes:
    @pytest.mark.parametrize("bad_genotype", ["--", "DD", "", "XX", "00", "ZZ", "NaN"])
    def test_malformed_does_not_crash(self, bad_genotype):
        """Pipeline should not raise for any malformed genotype value."""
        df = make_df([
            {"rsid": "rs_bad", "chrom": "3", "pos": 999, "genotype": bad_genotype},
            {"rsid": "rs_good", "chrom": "3", "pos": 1000, "genotype": "AG"},
        ])

        with patch("src.sanitize.load_blocklist", return_value={}), \
             patch("src.sanitize.parse_genome", return_value=df), \
             patch("src.sanitize.get_ld_neighbors", return_value=[]):

            result = sanitize("fake_path.txt")  # should not raise

        assert len(result) == 2

    def test_malformed_not_replaced(self):
        """A blocklisted SNP with '--' genotype should NOT get blocklist replacement."""
        df = make_df([
            {"rsid": "rs123", "chrom": "1", "pos": 100, "genotype": "--"},
        ])

        with patch("src.sanitize.load_blocklist", return_value=FAKE_BLOCKLIST), \
             patch("src.sanitize.parse_genome", return_value=df), \
             patch("src.sanitize.get_ld_neighbors", return_value=[]):

            result = sanitize("fake_path.txt")

        # Should keep "--", not replace with "AA" from blocklist
        assert result.iloc[0]["genotype"] == "--"

    def test_malformed_not_ld_tagged(self):
        """A malformed SNP that happens to be an LD neighbor should NOT be tagged."""
        df = make_df([
            {"rsid": "rs123", "chrom": "1", "pos": 100, "genotype": "CT"},
            {"rsid": "rs789", "chrom": "1", "pos": 300, "genotype": "--"},
        ])

        with patch("src.sanitize.load_blocklist", return_value={"rs123": "AA"}), \
             patch("src.sanitize.parse_genome", return_value=df), \
             patch("src.sanitize.get_ld_neighbors", return_value=["rs789"]):

            result = sanitize("fake_path.txt")

        rs789 = result[result["rsid"] == "rs789"].iloc[0]
        assert rs789["ld_neighbor"] == False
        assert rs789["genotype"] == "--"


# ------------------------------------------------------------------
# validate_and_clean unit tests (no mocking needed)
# ------------------------------------------------------------------

class TestValidateAndClean:
    def test_double_dash_flagged_as_skip(self):
        df = make_df([{"rsid": "rs1", "chrom": "1", "pos": 1, "genotype": "--"}])
        result = validate_and_clean(df)
        assert result.iloc[0]["_skip"] == True

    def test_empty_string_flagged_as_skip(self):
        df = make_df([{"rsid": "rs1", "chrom": "1", "pos": 1, "genotype": ""}])
        result = validate_and_clean(df)
        assert result.iloc[0]["_skip"] == True

    def test_valid_genotype_not_skipped(self):
        df = make_df([{"rsid": "rs1", "chrom": "1", "pos": 1, "genotype": "AG"}])
        result = validate_and_clean(df)
        assert result.iloc[0]["_skip"] == False

    def test_genotype_normalized_to_uppercase(self):
        df = make_df([{"rsid": "rs1", "chrom": "1", "pos": 1, "genotype": "ag"}])
        result = validate_and_clean(df)
        assert result.iloc[0]["genotype"] == "AG"