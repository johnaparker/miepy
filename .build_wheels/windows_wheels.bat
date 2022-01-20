:: Activate and install software
::     conda create -n py37 python=3.7
::     conda activate py37
::     conda install -c anaconda python pip numpy scipy matplotlib tqdm sympy pandas pyyaml cmake make
::     conda install -c conda-forge quaternion spherical gsl eigen
::     conda install -c vpython vpython

:: Manually build cpp with CMake
::     cmake -G "Visual Studio 15 2017 Win64" ..
::     cmake .. -A x64 --config Release
::     cmake --build . -- /m

set root=C:\Users\john7\Miniconda3
call %root%\Scripts\activate.bat %root%

:: call conda update --yes conda
:: call conda update --yes anaconda

for %%x in (35,36,37) do (
    call activate py%%x
    :: call conda update anaconda
    call conda install --yes -c anaconda pip numpy scipy matplotlib tqdm sympy pandas pyyaml cmake make
    call conda install --yes -c conda-forge quaternion spherical gsl eigen
    :: call conda update --all --yes
    call copy /y %root%\envs\py%%x\Library\bin\gsl.dll dlls
    call pip wheel . -w .build_wheels\wheelhouse --no-deps -v
)
