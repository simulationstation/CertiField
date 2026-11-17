#!/usr/bin/env wolframscript
(* Parse every .wl file of the paclet with CodeParser, without evaluating anything, and fail on syntax errors.
   Usage:  wolframscript -file tools/parse_check.wls            (from the repository root)          *)

Needs["CodeParser`"];

root = If[Length[$ScriptCommandLine] > 1, $ScriptCommandLine[[2]], Directory[]];
files = Join[FileNames["*.wl", FileNameJoin[{root, "Kernel"}]], FileNames["*.wl", root], FileNames["*.wlt", FileNameJoin[{root, "Tests"}]]];

problems[file_] := Module[{cst = CodeParser`CodeParse[File[file]], issues, nodes},
  issues = Cases[cst, (CodeParser`SyntaxIssue | CodeParser`EncodingIssue)[tag_, msg_, sev_, ___] /; MemberQ[{"Error", "Fatal"}, sev] :> tag <> ": " <> msg, Infinity];
  nodes = Cases[cst, (CodeParser`ErrorNode | CodeParser`UnterminatedGroupNode | CodeParser`UnterminatedCallNode | CodeParser`GroupMissingCloserNode)[tag_, ___] :> ToString[tag], Infinity];
  Join[issues, nodes]];

(* self-test: an unterminated nested comment must be detected *)
selfTest = With[{tmp = FileNameJoin[{$TemporaryDirectory, "certifield_parse_selftest.wl"}]},
  WriteString[tmp, "(* text with a (*) inside *)\nx = 1;\n"]; Close[tmp];
  With[{r = problems[tmp]}, DeleteFile[tmp]; r =!= {}]];
If[! selfTest, Print["parse_check: self-test FAILED (nested comments are not detected)"]; Exit[2]];

bad = Select[Table[{FileNameTake[f], problems[f]}, {f, files}], #[[2]] =!= {} &];
If[bad === {},
  Print["parse_check: OK, ", Length[files], " files parsed cleanly"]; Exit[0],
  Scan[Print["parse_check: ", #[[1]], ": ", Short[#[[2]], 3]] &, bad]; Exit[1]];
