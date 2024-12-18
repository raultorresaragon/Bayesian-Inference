# year dummies
dummify <- function(y, ad) {
  ad$temp <- 0
  ad$temp[ad$year == y] <- 1
  names(ad)[names(ad) == "temp"] <- paste0("y", y)
  ad
}
for(y in unique(ad$year)) {
  ad <- dummify(y, ad)
}

# state dummies
dummify <- function(s, ad) {
  ad$temp <- 0
  ad$temp[ad[["state"]] == s] <- 1
  names(ad)[names(ad) == "temp"] <- paste0("state_", s)
  ad
}
for(s in unique(ad$state)) {
  ad <- dummify(s, ad)
}