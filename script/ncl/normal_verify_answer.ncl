if (.not. isvar("fd1")) then
    print("Usage: ncl analyzefile.ncl " + "'fd1=/path/to/your/file' " + "'fd2=/path/to/your/file' ")
    exit
end if

f1 = addfile(fd1, "r")

f2 = addfile(fd2, "r")


; strange here
; if (isStrSubset(fd1, "adv")) then


    printMinMax((f1->tm - f2->tm), 0)
    print(sum(f1->tm - f2->tm) / sum(f1->tm))

    printMinMax((f1->te - f2->te), 0)
    print(sum(f1->te - f2->te) / sum(f1->te))

    printMinMax((f1->zs - f2->zs), 0)
    print(sum(f1->zs - f2->zs) / sum(f1->zs))

    printMinMax((f1->dzsdx - f2->dzsdx), 0)
    print(sum(f1->dzsdx - f2->dzsdx) / sum(f1->dzsdx))

    printMinMax((f1->dzsdy - f2->dzsdy), 0)
    print(sum(f1->dzsdy - f2->dzsdy) / sum(f1->dzsdy))

    printMinMax((f1->phs - f2->phs), 0)
    print(sum(f1->phs - f2->phs) / sum(f1->phs))


    printMinMax((f1->u - f2->u), 0)
    print(sum(f1->u - f2->u) / sum(f1->u))

    printMinMax((f1->v - f2->v), 0)
    print(sum(f1->v - f2->v) / sum(f1->v))
    
    printMinMax((f1->pt - f2->pt), 0)
    print(sum(f1->pt - f2->pt) / sum(f1->pt))
    
    printMinMax((f1->t - f2->t), 0)
    print(sum(f1->t - f2->t) / sum(f1->t))
   
    printMinMax((f1->z - f2->z), 0)
    print(sum(f1->z - f2->z) / sum(f1->z))
   
    printMinMax((f1->ph - f2->ph), 0)
    print(sum(f1->ph - f2->ph) / sum(f1->ph))
   
    printMinMax((f1->vor - f2->vor), 0)
    print(sum(f1->vor - f2->vor) / sum(f1->vor))
   
    printMinMax((f1->div - f2->div), 0)
    print(sum(f1->div - f2->div) / sum(f1->div))
   
    
; else if (isStrSubset(fd1, "swm")) then

; else
 
; end if