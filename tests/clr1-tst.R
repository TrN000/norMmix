require(norMmix)

## "large" "lambda" values :
lam <- 709 - (0:4)/10
(w <- clr1inv(lam))
stopifnot(length(w) == 1 + length(lam),
          w[1] == 0, w[-1] > 0,
          all.equal(sum(w), 1))
