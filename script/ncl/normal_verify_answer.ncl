if (.not. isvar("fd1")) then
    print("Usage: ncl analyzefile.ncl " + "'fd1=/path/to/your/file' " + "'fd2=/path/to/your/file' ")
    exit
end if

f1 = addfile(fd1, "r")

f2 = addfile(fd2, "r")


; strange here
; if (isStrSubset(fd1, "adv")) then


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

    print(sum(abs(f1->tm - f2->tm)) / sum(abs(f1->tm)))
    print(sum(abs(f1->te - f2->te)) / sum(abs(f1->te)))
    print(sum(abs(f1->zs - f2->zs)) / sum(abs(f1->zs)))
    print(sum(abs(f1->dzsdx - f2->dzsdx)) / sum(abs(f1->dzsdx)))
    print(sum(abs(f1->dzsdy - f2->dzsdy)) / sum(abs(f1->dzsdy)))
    print(sum(abs(f1->phs - f2->phs)) / sum(abs(f1->phs)))
    print(sum(abs(f1->u - f2->u)) / sum(abs(f1->u)))
    print(sum(abs(f1->v - f2->v)) / sum(abs(f1->v)))
    print(sum(abs(f1->pt - f2->pt)) / sum(abs(f1->pt)))
    print(sum(abs(f1->t - f2->t)) / sum(abs(f1->t)))
    print(sum(abs(f1->z - f2->z)) / sum(abs(f1->z)))
    print(sum(abs(f1->ph - f2->ph)) / sum(abs(f1->ph)))
    print(sum(abs(f1->vor - f2->vor)) / sum(abs(f1->vor)))
    print(sum(abs(f1->div - f2->div)) / sum(abs(f1->div)))

; else if (isStrSubset(fd1, "swm")) then

; else
 
; end if