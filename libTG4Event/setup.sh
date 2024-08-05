root -l -q -e 'std::cout<<Form("rootcint -v1 -f libTG4EventProjectDict.cxx %s libTG4EventProjectHeaders.h libTG4EventLinkDef.h", gSystem->GetIncludePath())<<endl;'>MAKEP

rootCode='
TString cmd = gSystem->GetMakeSharedLib();
cmd.ReplaceAll("$BuildDir",".");
cmd.ReplaceAll("$Opt",gSystem->GetFlagsOpt());
cmd.ReplaceAll("$IncludePath",Form("%s %s",gSystem->GetIncludePath(), "-IlibTG4Event"));
cmd.ReplaceAll("$SourceFiles","libTG4EventProjectSource.cxx");
cmd.ReplaceAll("$LinkedLibs",gSystem->GetLibraries("","SDL"));
cmd.ReplaceAll("$SharedLib","libTG4Event.so");
cmd.ReplaceAll("$ObjectFiles","libTG4EventProjectSource.o");
std::cout<<cmd;
'
echo "$rootCode" | root -l -b  >> MAKEP
source MAKEP
