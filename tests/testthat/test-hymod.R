test_that('hymod', {

  data("BassRiver")

  tmp = rgn(simFunc=simFunc_hymod,
            par=c(400.,0.5,0.1,0.2,0.1),
            lower=c(1.,0.1,0.05,0.000001,0.000001),
            upper=c(1000.,2.,0.95,0.99999,0.99999),
            simTarget=BassRiverData$Runoff.mm.day[365:length(BassRiverData$Date)],
            nWarmUp=365,
            rain=BassRiverData$Rain.mm,
            pet=BassRiverData$ET.mm,
            control=list(dump=0))
  error=tmp$error;message=tmp$message;par=tmp$par;value=tmp$value;counts=tmp$counts;info=tmp$info

  # expected output based on F90 RGN hymod example
  expect_equal(signif(par,digits=7),
               c(146.7564,0.3635988,0.1895957,0.99999,0.7430698))
  expect_equal(signif(value,digits=7),6840.165)
  expect_equal(counts,367)
  expect_equal(info$nIter,29)
  expect_equal(info$termFlag,2)

})
