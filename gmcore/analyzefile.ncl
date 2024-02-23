; fd="./adv_mv.360x180.h0.nc"
  if (.not. isvar("fd")) then
      print("Usage: ncl analyzefile.ncl " + "'fd=/path/to/your/file' ")
      exit
  end if
f1 = addfile(fd, "r")

; retrieve the global attributes from the file
fatt_names = getvaratts(f1)   
if(.not.all(ismissing(fatt_names))) then
    print("")
    print("The file's global attributes are: ")
    print("")
    
    do i = 0,dimsizes(fatt_names)-1
    print(fatt_names(i)+" = "+f1@$fatt_names(i)$)
    end do

end if

print("Variables in this file are:")

    vNames = getfilevarnames (f1) ; get names of all variables on file
        
    print (vNames)               ; print all variable names on file