#!/usr/bin/env wolframscript
(* Build CertiField-<version>.paclet from a clean staging copy (only the files a user needs), then install it into a
   throw-away paclet directory and smoke-test it.  Usage, from the repository root:
     wolframscript -file tools/build_paclet.wls                                                        *)

root = DirectoryName[$InputFileName, 2];
build = FileNameJoin[{root, "build"}];
stage = FileNameJoin[{build, "CertiField"}];
Quiet[DeleteDirectory[stage, DeleteContents -> True]];
CreateDirectory[stage, CreateIntermediateDirectories -> True];
Scan[CopyFile[FileNameJoin[{root, #}], FileNameJoin[{stage, #}]] &, {"PacletInfo.wl", "LICENSE", "README.md"}];
(* docs: CertiFieldModelInfo and the certificates refer to docs/SCALAR_ND.md *)
Scan[CopyDirectory[FileNameJoin[{root, #}], FileNameJoin[{stage, #}]] &, {"Kernel", "Data", "Examples", "docs"}];

archive = CreatePacletArchive[stage, build];
Print["built: ", archive, "  (", Round[FileByteCount[archive]/1024.], " KB)"];

(* smoke test in this kernel: install from the archive and use it *)
PacletInstall[archive, ForceVersionInstall -> True];
Needs["CertiField`"];
Print["installed from: ", CertiField`Private`$CertiFieldDirectory];
Print[CertifiedVariationalBound["Phi4Kink"]["Summary"]];
Print[CertifiedVariationalBound[CertiFieldModel["ScalarFieldND", <|"Potential" -> Function[u, u^2 - u], "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 16|>]]["Summary"]];
Print[CertifiedVariationalBound[CertiFieldModel["LinearElasticity", <|"Lame" -> {3/2, 1}, "Domain" -> {{0, 1}, {0, 1}}, "Cells" -> 8,
    "TractionFreeFaces" -> {{2, "Lower"}, {2, "Upper"}}, "BoundaryData" -> Function[{x, y}, {x/10, 0}]|>]]["Summary"]];
(* leave no installed copy behind (it would shadow a development checkout) *)
PacletUninstall["CertiField"];
