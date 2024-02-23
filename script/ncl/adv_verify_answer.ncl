if (.not. isvar("fd1")) then
    print("Usage: ncl analyzefile.ncl " + "'fd1=/path/to/your/file' " + "'fd2=/path/to/your/file' ")
    exit
end if

f1 = addfile(fd1, "r")

f2 = addfile(fd2, "r")


; strange here
; if (isStrSubset(fd1, "adv")) then
    printMinMax(f1->lat - f2->lat, 0)
    printMinMax(f1->lon - f2->lon, 0)
    printMinMax(f1->lev - f2->lev, 0)
    printMinMax(f1->ilon - f2->ilon, 0)
    printMinMax(f1->ilat - f2->ilat, 0)
    printMinMax(f1->ilev - f2->ilev, 0)
    printMinMax(f1->z - f2->z, 0)
    printMinMax(f1->q0 - f2->q0, 0)
    printMinMax(f1->q1 - f2->q1, 0)
; else if (isStrSubset(fd1, "swm")) then

; else
 
; end if