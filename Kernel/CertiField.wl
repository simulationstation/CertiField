(* ::Package:: *)

(* ::Title:: *)
(*CertiField: certified bounds for variational problems*)

(* ::Text:: *)
(*CertiField returns rigorous two-sided bounds on the minimal energy of supported nonlinear field theories,*)
(*over the whole continuous admissible class (not a discretisation of it).  A candidate from any solver gives the*)
(*upper bound; a proof-carrying lower-bound argument (calibration, topology, level sets) gives the lower bound;*)
(*interval and exact rational arithmetic make every number in the report rigorous.  When a bound cannot be*)
(*established the result says "Inconclusive" instead of reporting an uncertified number.*)

BeginPackage["CertiField`"];

CertifiedVariationalBound::usage =
  "CertifiedVariationalBound[model, candidate] returns a CertifiedBound object with rigorous bounds \
LowerBound <= inf E <= UpperBound for the minimal energy of model over its admissible class, the upper bound \
coming from candidate.\nCertifiedVariationalBound[model] uses the model's built-in reference candidate.\n\
Options: \"GapTolerance\", \"Verbose\"; \"Resolution\" (SkyrmeDegreeOne); \"Nodes\", WorkingPrecision (one-dimensional \
models); Method (ScalarFieldND, SkyrmeDegreeOne), \"RefineCandidate\" (ScalarFieldND, LinearElasticity). CertiFieldModelInfo[name] \
lists the parameters of each model.";

CertifiedBound::usage =
  "CertifiedBound[...] is the result of CertifiedVariationalBound. Query it with cert[\"LowerBound\"], \
cert[\"UpperBound\"], cert[\"RelativeGap\"], cert[\"Status\"], cert[\"Theorem\"], cert[\"Properties\"].";

CertiFieldModel::usage =
  "CertiFieldModel[name, params] specifies a supported variational model. Supported names are listed in \
$CertiFieldModels; CertiFieldModelInfo[name] describes each one.";

VerifyCertificate::usage =
  "VerifyCertificate[cert] independently re-runs every rigorous step recorded in the CertifiedBound cert and \
returns True if all bounds are reproduced.";

CertiFieldModelInfo::usage =
  "CertiFieldModelInfo[name] returns an Association describing the energy, admissible class, theorem and \
parameters of a supported model.";

$CertiFieldModels::usage = "$CertiFieldModels lists the model names supported by CertifiedVariationalBound.";
$CertiFieldVersion::usage = "$CertiFieldVersion is the version of the installed CertiField paclet.";

CertifiedVariationalBound::nomodel = "`1` is not a supported model. Supported models: `2`.";
CertifiedVariationalBound::badparam = "`1`";
CertifiedVariationalBound::badcand = "`1`";

Begin["`Private`"];

$CertiFieldVersion = "0.3.0";
$CertiFieldDirectory = DirectoryName[$InputFileName, 2];

failure[tag_String, msg_String, extra_Association : <||>] :=
  Failure[tag, Join[<|"MessageTemplate" -> msg|>, extra]];

Get[FileNameJoin[{DirectoryName[$InputFileName], #}]] & /@ {
  "RigorousArithmetic.wl",
  "RigorousQuadrature.wl",
  "Minorants.wl",
  "SkyrmeLevelSet.wl",
  "SkyrmeHedgehog.wl",
  "SkyrmeAnyMass.wl",
  "SkyrmeLevelSetDegree.wl",
  "ScalarField1D.wl",
  "ScalarFieldND.wl",
  "ThinFilmMagnet.wl",
  "Elasticity.wl",
  "CertifiedBound.wl",
  "Models.wl"};

End[];
EndPackage[];
