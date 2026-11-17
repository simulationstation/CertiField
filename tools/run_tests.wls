#!/usr/bin/env wolframscript
(* Run the CertiField test suites from the repository root.
     wolframscript -file tools/run_tests.wls            fast suite
     wolframscript -file tools/run_tests.wls --skyrme   fast suite + Skyrme suite (a few minutes)   *)

root = DirectoryName[$InputFileName, 2];
PacletDirectoryLoad[root];
files = {FileNameJoin[{root, "Tests", "CertiField.wlt"}]};
If[MemberQ[$ScriptCommandLine, "--skyrme"], AppendTo[files, FileNameJoin[{root, "Tests", "Skyrme.wlt"}]]];

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
