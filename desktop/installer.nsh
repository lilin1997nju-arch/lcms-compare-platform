!macro customRemoveFiles
  IfFileExists "$INSTDIR.lcms-data-preserved" 0 +2
    Abort "A preserved LC-MS data directory already exists beside the installation directory. Please recover or rename it before continuing."

  ClearErrors
  Rename "$INSTDIR\data" "$INSTDIR.lcms-data-preserved"
  ClearErrors
  Rename "$INSTDIR\lcms-desktop-config.json" "$INSTDIR.lcms-config-preserved"

  RMDir /r "$INSTDIR"
  CreateDirectory "$INSTDIR"

  ClearErrors
  Rename "$INSTDIR.lcms-data-preserved" "$INSTDIR\data"
  ClearErrors
  Rename "$INSTDIR.lcms-config-preserved" "$INSTDIR\lcms-desktop-config.json"
!macroend
