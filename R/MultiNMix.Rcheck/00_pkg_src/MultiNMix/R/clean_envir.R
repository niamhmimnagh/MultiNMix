#' @title Clean Environment of Temporary Objects
#'
#' @description This function removes temporary objects from a specified environment. It is designed for use in contexts where intermediate objects are created and should be cleaned up after execution to prevent clutter or unintended behavior.
#'
#' @details This is an auxiliary function that removes temporary objects created by Nimble from the environment.
#'
#' @param env_rm An environment from which to remove objects. Defaults to the global environment (`.GlobalEnv`). This may be set to a different environment for more controlled cleanup.
#'
#' @return The function does not return any values.
#'
#' @export
#'
clean_envir <- function(env_rm = .GlobalEnv) {
  # # Debugging output to check contents before cleanup
  # cat("Before cleanup:", ls(env_rm), "\n")

  # List all objects that start with "str" and confirm they exist
  str_objects <- ls(envir = env_rm, pattern = "^str")

  # If any objects match, remove them
  if (length(str_objects) > 0) {
    rm(list = str_objects, envir = env_rm)
  }

  # # Debugging output to confirm cleanup
  # cat("After cleanup:", ls(env_rm), "\n")
}
