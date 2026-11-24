#!/usr/bin/env wolframscript
(* Run the CertiField test suites from the repository root.
     wolframscript -file tools/run_tests.wls            fast suite
     wolframscript -file tools/run_tests.wls --skyrme   fast suite + Skyrme suite (a few minutes)
     wolframscript -file tools/run_tests.wls --skyrme-any-mass   fast suite + Skyrme for a mass that is not shipped (about eight minutes)
     wolframscript -file tools/run_tests.wls --nd       fast suite + ScalarFieldND suite (about two minutes)
     wolframscript -file tools/run_tests.wls --nd-slow  fast suite + the long ScalarFieldND runs (about ten minutes)
     wolframscript -file tools/run_tests.wls --magnet   fast suite + ThinFilmMagnet suite (about two minutes)
     wolframscript -file tools/run_tests.wls --elasticity  fast suite + LinearElasticity suite (under a minute)
     wolframscript -file tools/run_tests.wls --skyrme-degree   fast suite + Skyrme fields of degree two (about six minutes)   *)

root = DirectoryName[$InputFileName, 2];
PacletDirectoryLoad[root];
files = {FileNameJoin[{root, "Tests", "CertiField.wlt"}]};
If[MemberQ[$ScriptCommandLine, "--skyrme"], AppendTo[files, FileNameJoin[{root, "Tests", "Skyrme.wlt"}]]];
If[MemberQ[$ScriptCommandLine, "--skyrme-any-mass"], AppendTo[files, FileNameJoin[{root, "Tests", "SkyrmeAnyMass.wlt"}]]];
If[MemberQ[$ScriptCommandLine, "--nd"], AppendTo[files, FileNameJoin[{root, "Tests", "ScalarFieldND.wlt"}]]];
If[MemberQ[$ScriptCommandLine, "--nd-slow"], AppendTo[files, FileNameJoin[{root, "Tests", "ScalarFieldNDSlow.wlt"}]]];
If[MemberQ[$ScriptCommandLine, "--magnet"], AppendTo[files, FileNameJoin[{root, "Tests", "ThinFilmMagnet.wlt"}]]];
If[MemberQ[$ScriptCommandLine, "--elasticity"], AppendTo[files, FileNameJoin[{root, "Tests", "Elasticity.wlt"}]]];
If[MemberQ[$ScriptCommandLine, "--skyrme-degree"], AppendTo[files, FileNameJoin[{root, "Tests", "SkyrmeDegree.wlt"}]]];

failed = 0;
Do[
  Module[{report = TestReport[f]},
    Print[FileNameTake[f], ": ", report["TestsSucceededCount"], " passed, ", report["TestsFailedCount"], " failed  (",
      Round[QuantityMagnitude[report["TimeElapsed"]], 0.1], " s)"];
    Scan[Print["  FAILED ", #["TestID"], ":  ", Short[#["ActualOutput"], 2]] &,
      Join @@ Values[report["TestsFailed"]]];
    failed += report["TestsFailedCount"]],
  {f, files}];
Exit[If[failed == 0, 0, 1]];
