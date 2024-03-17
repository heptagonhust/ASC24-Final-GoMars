if (.not. isvar("fd1")) then
    print("Usage: ncl analyzefile.ncl " + "'fd1=/path/to/your/file' " + "'fd2=/path/to/your/file' ")
    exit
end if

f1 = addfile(fd1, "r")

f2 = addfile(fd2, "r")


; strange here
; if (isStrSubset(fd1, "adv")) then

    printMinMax((f1->lat - f2->lat) / f2->lat, 0)

    averagelon = avg(f1->lon)
    printMinMax(f1->lon, 0)
    printMinMax(f1->lon - f2->lon, 0)
    print(averagelon)

    averagelev = avg(f1->lev)
    printMinMax(f1->lev, 0)
    printMinMax(f1->lev - f2->lev, 0)
    print(averagelev)

    averageilon = avg(f1->ilon)
    printMinMax(f1->ilon, 0)
    printMinMax(f1->ilon - f2->ilon, 0)
    print(averageilon)

    averageilat = avg(f1->ilat)
    printMinMax(f1->ilat - f2->ilat, 0)
    print(averageilat)

    printMinMax(f1->ilev - f2->ilev, 0)
    averageilev = avg(f1->ilev)
    print(averageilev)

    printMinMax(f1->tm - f2->tm, 0)
    averagetm = avg(f1->tm)
    print(averagetm)

    printMinMax(f1->te - f2->te, 0)
    averagete = avg(f1->te)
    print(averagete)

    printMinMax(f1->zs - f2->zs, 0)
    averagezs = avg(f1->zs)
    print(averagezs)

    printMinMax(f1->dzsdx - f2->dzsdx, 0)
    averagedzsdx = avg(f1->dzsdx)
    print(averagedzsdx)

    printMinMax(f1->dzsdy - f2->dzsdy, 0)
    averagedzsdy = avg(f1->dzsdy)
    print(averagedzsdy)

    printMinMax(f1->phs - f2->phs, 0)
    averagephs = avg(f1->phs)
    print(averagephs)

    printMinMax(f1->u - f2->u, 0)
    averageu = avg(f1->u)
    print(averageu)

    printMinMax(f1->v - f2->v, 0)
    averagev = avg(f1->v)
    print(averagev)
    
    printMinMax(f1->pt - f2->pt, 0)
    averagept = avg(f1->pt)
    print(averagept)
    
    printMinMax(f1->t - f2->t, 0)
    averaget = avg(f1->t)
    print(averaget)
   
    printMinMax(f1->z - f2->z, 0)
    averagez = avg(f1->z)
    print(averagez)
   
    printMinMax(f1->ph - f2->ph, 0)
    averageph = avg(f1->ph)
    print(averageph)
   
    printMinMax(f1->vor - f2->vor, 0)
    averagevor = avg(f1->vor)
    print(averagevor)
   
    printMinMax(f1->div - f2->div, 0)
    averagediv = avg(f1->div)
    print(averagediv)
   
    printMinMax(f1->landmask - f2->landmask, 0)
    averageland = avg(f1->landmask)
    print(averageland)

; else if (isStrSubset(fd1, "swm")) then

; else
 
; end if