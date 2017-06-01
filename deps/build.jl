if is_unix()
    
    try
       run(`sudo pip install gcvspline`)
    catch
       error("pip not install. Please install pip")
    end
else # windows
    try
       run(`pip install gcvspline`)
    catch
       error("pip not install. Please install pip")
    end
end
	
# Windows example from Dierckx.jl... Don't have windows so I can't really work on this for GCVSPL

#@windows_only begin
#    # these binaries were cross-compiled from Cygwin via the following steps:
#    # mkdir -p bin32 && mkdir -p bin64
#    # i686-w64-mingw32-gfortran -o bin32/libddierckx.dll -O3 -shared \
#    #   -static-libgfortran -static-libgcc src/ddierckx/*.f
#    # x86_64-w64-mingw32-gfortran -o bin64/libddierckx.dll -O3 -shared \
#    #   -static-libgfortran -static-libgcc src/ddierckx/*.f
#    url = "https://cache.julialang.org/https://bintray.com/artifact/download/tkelman/generic/ddierckx.7z"
#    try
#        run(`curl -LO $url`)
#    catch
#        run(`powershell -Command "(new-object net.webclient).DownloadFile(\"$url\", \"ddierckx.7z\")"`)
#    end
#    run(`7z x -y ddierckx.7z`)
#end