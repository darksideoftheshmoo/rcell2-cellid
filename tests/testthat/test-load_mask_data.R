test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("read_coords_tsv works on examples", {
  # Example: Choose one set of example images:
  data.dir <- system.file("extdata/sample_datasets/sample_time_series/",
                          package = "rcell2.examples")
  # Load X
  # X <- cell.load.alt(path = data.dir)
  # saveRDS(X, "inst/extdata/test_data/cell.load.alt.RDS")  # Save X ref
  
  # Load X ref
  ref.X <- system.file("extdata/test_data/cell.load.alt.RDS",
                       package = "rcell2.cellid")
  X <- readRDS(ref.X)
  
  # Load Y
  Y <- cell.load.boundaries(data.source = "masks.tsv", positions = X$positions)
  # saveRDS(Y, "inst/extdata/test_data/cell.load.boundaries.RDS")  # Save Y ref
  
  # Load Y ref
  Y.ref <- system.file("extdata/test_data/cell.load.boundaries.RDS",
                         package = "rcell2.cellid")
  Y.ref <- readRDS(Y.ref)
  
  # Test
  expect_equal(Y, Y.ref)
})
