@echo OFF
set gbfilter="D:\Documents\Projets\Isotropix\BGFilter\Debug\BGFilter.exe"
set inputPath="D:\Documents\Projets\Isotropix\BGFilter\"
set Start
set End
echo Test Starting Time: %TIME%
for /l %%H in (1, 9, 100) do (
	for /l %%W in (1, 9, 100) do (
		for /l %%K in (0, 1, 5) do (
			for /l %%L in (0, 1, 9) do (
				for %%I in (*_in.bmp) do (
					echo %%H %%W %%K.%%L %%I
					echo   Start: %TIME%
					%gbfilter% "%%~fI" "%%~pI%%~nI_%%H_%%W_%%K.%%L_out%%~xI" %%K.%%L %%W %%H
					echo   End  : %TIME%
					REM set /A StartTime=%Start:~0,2%*3600000 + %Start:~3,2%*60000 + %Start:~6,2%*1000 + %Start:~9,2%*10
					REM set /A EndTime=%End:~0,2%*3600000 + %End:~3,2%*60000 + %End:~6,2%*1000 + %End:~9,2%*10
					REM echo %StartTime% %EndTime%
					REM set /A Elapsed=%EndTime%-%StartTime%
				)
			)
		)
	)
)
echo Test Ending Time  : %TIME%
pause