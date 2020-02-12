"""

    bellerophon

Given exact integers f and e with f > 0,
return the float with n bits of precision that best approximates it.
Tries to do the calculation using floats with p bits of precision.
The error bounds used here assume perfect floating point arithmetic,
as in the IEEE standard. They are independent of p and n.
"""
function bellerophon(f::Integer, e::Integer; n::Int=53, p::Int=64, h::Int=2^10)
    log5of2n = Int(ceil(n / log2(5))) # = 23 for n = 53
    if f < 2^n && e >= 0 && e < h && e < log5of2n
        T(f) * exp10(e) # multiplication with n bits precision
    elseif f < 2^n && e < 0 && -e < h && -e < log5of2n
        T(f) / exp10(-e) # division with n bits of precision
    elseif f < 2^p && e >= 0 && e < h
        multiply_and_test(f, e, 0)
    elseif f < 2^p && (e < 0 || e >= h)
        multiply_and_test(f, e, 3)
    elseif f >= 2^p && e >= 0 && e < h
        multiply_and_test(f, e, 1)
    elseif f >= 2^p && (e < 0 || e >= h)
        multiply_and_test(f, e, 4)
    end
end
#=        
(define (Bellerophon f e)
  (cond ((and (< f two^n) (>= e 0) (< e h) (< e log5-of-two^n))
         (shortfloat-multiply (int->shortfloat f)
                              (float->shortfloat (ten-to-e e))))
        ((and (< f two^n) (< e 0) (< (- e) h) (< (- e) log5-of-two^n))
         (shortfloat-divide (int->shortfloat f)
                            (float->shortfloat (ten-to-e (- e)))))
        ((and (< f two^p) (>= e 0) (< e h))
         (multiply-and-test f e 0))
        ((and (< f two^p) (or (< e 0) (>= e h)))
         (multiply-and-test f e 3))
        ((and (>= f two^p) (>= e 0) (< e h))
         (multiply-and-test f e 1))
        ((and (>= f two^p) (or (< e 0) (>= e h)))
         (multiply-and-test f e 4))))
=#
"""
Slop, expressed in units of the least significant bit, is an
inclusive bound for the error accumulated during the floating
point calculation of an approximation to f * 10^e. (Slop is
not a bound for the true error, but bounds the difference
between the approximation z and the best possible approximation
that uses p bits of significand.)
Fail is a slow but perfect backup algorithm.
Z is passed so fail can use it as a starting approximation.
"""
define multiply_and_test(f::Integer, e::Int, slop::Int)
    x = float(f)
    y = exp10(e)
    z = mult(x, y)


#=
(define (multiply-and-test f e slop)
  (let ((x (int->float f))
        (y (ten-to-e e)))
    (let ((z (float-multiply x y)))
      (let ((lowbits (remainder (float-significand z) two^p-n)))

        ; is the slop large enough to make a difference when
        ; rounding to n bits?

        (if (<= (abs (- lowbits two^p-n-1)) slop)
            (fail f e z)
            (float->shortfloat z))))))

(define (fail f e z)
  (AlgorithmR f e (float->shortfloat z)))

(define n 53)   ; IEEE double
(define p 54)   ; an extended precision
(define two^p (expt 2 p))
(define two^p-1 (expt 2 (- p 1)))
(define two^p-n (expt 2 (- p n)))
(define two^p-n-1 (expt 2 (- p n 1)))
(define two^n (expt 2 n))
(define two^n-1 (expt 2 (- n 1)))

(define log5-of-two^n
  (inexact->exact (ceiling (/ (log two^n) (log 5)))))
=#

