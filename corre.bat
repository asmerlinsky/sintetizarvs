gcc sintetizadorprueba.c -lm -o sintetizar
if %errorlevel% neq 0 exit /b %errorlevel%

sintetizar BAdia.sound BAcanto.Sound
if %errorlevel% neq 0 exit /b %errorlevel%

sintetizar BAnoche.sound BAcanto.Sound
if %errorlevel% neq 0 exit /b %errorlevel%

py graficarsintetico.py BAdia.sound BAcanto.Sound BAnoche.sound
