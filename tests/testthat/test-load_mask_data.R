test_that("cell.load.alt works on examples", {
  # Example: Choose one set of example images:
  data.dir <- system.file("extdata/sample_datasets/sample_time_series/",
                          package = "rcell2.examples")
  # Load X
  X <- cell.load.alt(path = data.dir)
  # saveRDS(X, "inst/extdata/test_data/cell.load.alt.RDS")  # Save X ref
  
  # Load X ref
  X.ref <- system.file("extdata/test_data/cell.load.alt.RDS",
                       package = "rcell2.cellid")
  X.ref <- readRDS(ref.X)
  
  # Test
  expect_equal(X, X.ref)
})

test_that("cell.load.boundaries works on examples", {
  # Example: Choose one set of example images:
  data.dir <- system.file("extdata/sample_datasets/sample_time_series/",
                          package = "rcell2.examples")
  # Load X
  # X <- cell.load.alt(path = data.dir)
  # saveRDS(X, "inst/extdata/test_data/cell.load.alt.RDS")  # Save X ref
  
  # Load X ref
  X.ref <- system.file("extdata/test_data/cell.load.alt.RDS",
                       package = "rcell2.cellid")
  X.ref <- readRDS(X.ref)
  
  # Load Y
  Y <- cell.load.boundaries(data.source = "masks.tsv", positions = X.ref$positions)
  # saveRDS(Y, "inst/extdata/test_data/cell.load.boundaries.RDS")  # Save Y ref
  
  # Load Y ref
  Y.ref <- system.file("extdata/test_data/cell.load.boundaries.RDS",
                         package = "rcell2.cellid")
  Y.ref <- readRDS(Y.ref)
  
  # Test
  expect_equal(Y, Y.ref)
})
