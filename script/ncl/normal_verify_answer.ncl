if (.not. isvar("fd1")) then
    print("Usage: ncl analyzefile.ncl " + "'fd1=/path/to/your/file' " + "'fd2=/path/to/your/file' ")
    exit
end if

f1 = addfile(fd1, "r")

f2 = addfile(fd2, "r")


; strange here
; if (isStrSubset(fd1, "adv")) then
    printMinMax(abs(f1->lat - f2->lat), 0)
    printMinMax(abs(f1->lon - f2->lon), 0)
    printMinMax(abs(f1->lev - f2->lev), 0)
    printMinMax(abs(f1->ilon - f2->ilon), 0)
    printMinMax(abs(f1->ilat - f2->ilat), 0)
    printMinMax(abs(f1->ilev - f2->ilev), 0)
    printMinMax(abs(f1->tm - f2->tm), 0)
    printMinMax(abs(f1->te - f2->te), 0)
    printMinMax(abs(f1->zs - f2->zs), 0)
    printMinMax(abs(f1->dzsdx - f2->dzsdx), 0)
    printMinMax(abs(f1->dzsdy - f2->dzsdy), 0)
    printMinMax(abs(f1->phs - f2->phs), 0)
    printMinMax(abs(f1->u - f2->u), 0)
    printMinMax(abs(f1->v - f2->v), 0)
    printMinMax(abs(f1->pt - f2->pt), 0)
    printMinMax(abs(f1->t - f2->t), 0)
    printMinMax(abs(f1->z - f2->z), 0)
    printMinMax(abs(f1->ph - f2->ph), 0)
    printMinMax(abs(f1->vor - f2->vor), 0)
    printMinMax(abs(f1->div - f2->div), 0)
    printMinMax(abs(f1->landmask - f2->landmask), 0)
; else if (isStrSubset(fd1, "swm")) then

; else
 
; end if