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
    printMinMax(f1->ilon - f2->ilon, 0)
    printMinMax(f1->ilat - f2->ilat, 0)
    printMinMax(f1->tm - f2->tm, 0)
    printMinMax(f1->te - f2->te, 0)
    printMinMax(f1->tpe - f2->tpe, 0)
    printMinMax(f1->u - f2->u, 0)
    printMinMax(f1->v - f2->v, 0)
    printMinMax(f1->z - f2->z, 0)
    printMinMax(f1->pv - f2->pv, 0)
    printMinMax(f1->div - f2->div, 0)
    printMinMax(f1->zs - f2->zs, 0)
; else if (isStrSubset(fd1, "swm")) then

; else
 
; end if