#ifndef AppVersion
  #error AppVersion must be supplied by the build script.
#endif
#ifndef BundleDir
  #error BundleDir must be supplied by the build script.
#endif
#ifndef ArtifactDir
  #error ArtifactDir must be supplied by the build script.
#endif

[Setup]
AppId={{A5AF6196-8A2A-45F6-B6CA-63701E0EEC37}
AppName=PichAnalysis
AppVersion={#AppVersion}
VersionInfoVersion={#AppVersion}.0
VersionInfoProductName=PichAnalysis
VersionInfoDescription=PichAnalysis
DefaultDirName={localappdata}\Programs\PichAnalysis
DefaultGroupName=PichAnalysis
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir={#ArtifactDir}
OutputBaseFilename=PichAnalysis-{#AppVersion}-Windows-x86_64-Setup
Compression=lzma2
SolidCompression=yes
UninstallDisplayIcon={app}\PichAnalysis.exe
DisableProgramGroupPage=yes
WizardStyle=modern
LanguageDetectionMethod=none

[Languages]
Name: "en"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: unchecked

[Files]
Source: "{#BundleDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\PichAnalysis"; Filename: "{app}\PichAnalysis.exe"
Name: "{autodesktop}\PichAnalysis"; Filename: "{app}\PichAnalysis.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\PichAnalysis.exe"; Description: "Launch PichAnalysis"; Flags: nowait postinstall skipifsilent
