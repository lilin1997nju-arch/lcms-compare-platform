!macro customRemoveFiles
  IfFileExists "$INSTDIR.lcms-data-preserved" 0 lcms_check_config_backup
    Abort "A preserved LC-MS data directory already exists beside the installation directory. Please recover or rename it before continuing."
  lcms_check_config_backup:
  IfFileExists "$INSTDIR.lcms-config-preserved" 0 lcms_preserve_data
    Abort "A preserved LC-MS configuration already exists beside the installation directory. Please recover or rename it before continuing."

  lcms_preserve_data:
  IfFileExists "$INSTDIR\data" 0 lcms_preserve_config
  ClearErrors
  Rename "$INSTDIR\data" "$INSTDIR.lcms-data-preserved"
  IfErrors 0 lcms_preserve_config
    Abort "Unable to preserve LC-MS data. Close the application and retry. Installation files have not been removed."

  lcms_preserve_config:
  IfFileExists "$INSTDIR\lcms-desktop-config.json" 0 lcms_remove_program
  ClearErrors
  Rename "$INSTDIR\lcms-desktop-config.json" "$INSTDIR.lcms-config-preserved"
  IfErrors 0 lcms_remove_program
    Rename "$INSTDIR.lcms-data-preserved" "$INSTDIR\data"
    Abort "Unable to preserve LC-MS configuration. Installation files have not been removed. Any preserved data remains beside the installation directory."

  lcms_remove_program:
  ClearErrors
  RMDir /r "$INSTDIR"
  IfErrors 0 lcms_restore_data
    Abort "Unable to remove old program files. Preserved data/configuration remain beside the installation directory for recovery."
  lcms_restore_data:
  CreateDirectory "$INSTDIR"

  IfFileExists "$INSTDIR.lcms-data-preserved" 0 lcms_restore_config
  ClearErrors
  Rename "$INSTDIR.lcms-data-preserved" "$INSTDIR\data"
  IfErrors 0 lcms_restore_config
    Abort "Unable to restore LC-MS data. Data is preserved beside the installation directory; recover it before continuing."
  lcms_restore_config:
  IfFileExists "$INSTDIR.lcms-config-preserved" 0 lcms_restore_complete
  ClearErrors
  Rename "$INSTDIR.lcms-config-preserved" "$INSTDIR\lcms-desktop-config.json"
  IfErrors 0 lcms_restore_complete
    Abort "Unable to restore LC-MS configuration. It is preserved beside the installation directory."
  lcms_restore_complete:
!macroend
