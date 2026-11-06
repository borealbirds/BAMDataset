# Writes a data.frame to a parquet file using DuckDB (no arrow dependency needed)
#
# The file is first written to "<path>.tmp" and then renamed, so a crashed or interrupted job never
# leaves a partially written file at "path".
#
# df: data.frame to write
# path: output file path (should end in .parquet)
#
# Returns path (invisibly)
write_parquet_duckdb = function(df,
                                path) {

  conn = DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(conn, shutdown = TRUE))

  tmp_path = paste0(path, ".tmp")
  duckdb::duckdb_register(conn, "df_to_write", df)
  DBI::dbExecute(conn, paste0("COPY df_to_write TO ", DBI::dbQuoteString(conn, tmp_path), " (FORMAT parquet)"))
  duckdb::duckdb_unregister(conn, "df_to_write")

  if (file.exists(path)) file.remove(path)
  file.rename(tmp_path, path)

  invisible(path)

}
