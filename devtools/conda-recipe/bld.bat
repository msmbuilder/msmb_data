mkdir "%PREFIX%\share\"
xcopy msmb_data "%PREFIX%\share\msmb_data" /e
if errorlevel 1 exit 1
