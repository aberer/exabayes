#! /usr/bin/Rscript

## this is a failed attempt at computing the boundaries for a partial
## sliding window proposal

 
max=10
min = 0.1

a = c(1,2,3,4,5,.2,1)
arel = a / sum(a)





getDelta = function(rates, ind, indo)
    {
        inc = range(rates[ind])
        mininc = inc[1]
        maxinc = inc[2]

        inco = range(rates[indo])
        minex = inco[1]
        maxex = inco[2]
        
        rs = rates[length(rates)]
        N = length(rates)
        n = length(ind)

        rhoi = max
        ri = maxinc
        d1 = (n * ri - n * rhoi * rs) / (rhoi - 1 )

        rhoi = min
        ri = mininc
        d2 = (n * ri - n * rhoi * rs) / (rhoi - 1 )
        
        rhoo = min
        ri = minex
        d3 =  (( - rhoo * rs + ri ) * (N - n)) / ( (rhoo * (N-n)) - 1 )

        rhoo = max
        ri = maxex
        d4 =  (( - rhoo * rs + ri ) * (N - n)) / ( (rhoo * (N-n)) - 1 )
        
        abs(c(d1,d2,d3,d4) )
    }


## getDelta2 = function(rates, ind, indo)
##     {
##         inc = range(rates[ind])
##         mininc = inc[1]
##         maxinc = inc[2]

##         inco = range(rates[indo])
##         minex = inco[1]
##         maxex = inco[2]
        
##         rs = rates[length(rates)]
##         N = length(rates)
##         n = length(ind)

##         rhoi = max
##         ri = maxinc
##         d1 = (n * (N - n ) * (rhoi * rs - ri)) / (N - n + n * rhoi)
        
##         rhoi = min
##         ri = mininc
##         d2 = (n * (N - n ) * (rhoi * rs - ri)) / (N - n + n * rhoi)
        
##         rhoo = min
##         ri = minex
##         d3 = ( (N-n) * (ri - rhoo* rs )) / (1-rhoo)

##         rhoo = max
##         ri = maxex
##         d4 = ( (N-n) * (ri - rhoo* rs )) / (1-rhoo)

##         abs(c(d1,d2,d3,d4))
##     }




applyDelta = function (rates, delta ,ind, indo, dir)
{
    n = length(ind)
    N = length(rates)

    if(dir)
        {
            rates[ind] =  rates[ind] + delta / n
            rates[indo] = rates[indo] - delta / (N - n)
        }
    else
        {
            rates[ind] =  rates[ind] - delta / n
            rates[indo] = rates[indo] + delta / (N - n)
        }

    rates / rates[length(rates)]
}

res1 = min(getDelta(arel, 1:3, 4:7 ))

b = applyDelta(arel, res1, 1:3,4:7 , T)
c = applyDelta(arel, res1, 1:3,4:7 , F)
