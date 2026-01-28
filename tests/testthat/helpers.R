### Typical usage:
### 
###     expect_error2(some op , "not supported")
### 
### This won't break like expect_error(some op, "not supported") does if
### the error message happens to contain other white spaces like \n's or
### \t's instead of a single space between "not" and "supported". These
### can be introduced in an unpredictable way by 'stop(wmsg(...))'.
expect_error2 <- function(object, expected_string, ...)
{
    regexp <- gsub("[\\s]+", "[\\\\s]+", expected_string, perl=TRUE)
    expect_error(object, regexp=regexp, perl=TRUE, ...)
}

