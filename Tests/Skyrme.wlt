(* ::Package:: *)

(* Skyrme test suite (a few minutes: builds and certifies the minorant table, runs the Quick certificate,
   then re-verifies it).  Run with:  wolframscript -file tools/run_tests.wls --skyrme *)

Needs["CertiField`"];

sk = CertifiedVariationalBound["SkyrmeDegreeOne", "PionMass" -> 1];

VerificationTest[sk["Status"], "Certified", TestID -> "skyrme-mu1-status"]

(* the independent Python verifier (levelset_cert3.py) gives 1.414023462272984 on the same grid *)
VerificationTest[Abs[sk["LowerBound"] - 1.414023462272984] < 10^-9, True, TestID -> "skyrme-mu1-matches-python-verifier"]

(* numerically converged hedgehog: E/(12 Pi^2) = 1.41622 at mu = 1 *)
VerificationTest[1.4162 < sk["UpperBound"] < 1.41623, True, TestID -> "skyrme-mu1-hedgehog-upper-bound"]
VerificationTest[sk["RelativeGap"] < 0.0016, True, TestID -> "skyrme-mu1-gap"]
VerificationTest[sk["Certificate"]["LevelSet"]["Pairs"], 5025904, TestID -> "skyrme-mu1-same-transitions-as-python"]
VerificationTest[VerifyCertificate[sk], True, TestID -> "skyrme-mu1-verifies"]
