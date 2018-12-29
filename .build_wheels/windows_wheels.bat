:: Activate and install software
::     conda activate py37
::     conda install -c anaconda cmake make
::     conda install -c conda-forge gsl eigen

:: Manually build cpp with CMake
::     cmake -G "Visual Studio 15 2017 Win64" ..
::     cmake .. -A x64 --config Release
::     cmake --build . -- /m

set root=C:\Users\john7\Miniconda3
call %root%\Scripts\activate.bat %root%

call conda update --yes -n root conda

for %%x in (35,36,37) do (
    call activate py%%x
    :: call conda update anaconda
    call conda install --yes -c anaconda cmake make
    call conda install --yes -c conda-forge gsl eigen
    call pip wheel . -w .build_wheels/wheelhouse --no-deps
)