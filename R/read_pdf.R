# Load required package (commented out to prevent execution during package load)
# if (!require(pdftools)) install.packages("pdftools")
# library(pdftools)

# 🔧 SET YOUR FOLDER PATH HERE (where your PDF books are stored)
# pdf_folder <- "C:/Users/tralucck/OneDrive"  # ← CHANGE THIS

# Create an output folder for text files
# output_folder <- file.path(pdf_folder, "extracted_texts")
# dir.create(output_folder, showWarnings = FALSE)

# Get all PDF files (case-insensitive)

# pdf_files <- list.files(
#   path = pdf_folder,
#   pattern = "\\.pdf$",
#   ignore.case = TRUE,
#   full.names = TRUE
# )

# Check if any PDFs were found
# if (length(pdf_files) == 0) {
#   stop("No PDF files found in: ", pdf_folder)
# }

# cat("Found", length(pdf_files), "PDF file(s). Starting extraction...\n\n")

# 🔁 Loop through each PDF (commented out to prevent execution during package load)
# for (pdf_path in pdf_files) {
#   filename <- basename(pdf_path)
#   cat("📄 Processing:", filename, "\n")
# 
#   # Generate output .txt path
#   txt_filename <- tools::file_path_sans_ext(filename)  # removes .pdf
#   txt_path <- file.path(output_folder, paste0(txt_filename, ".txt"))
# 
#   # Extract text with error handling
#   tryCatch({
#     text_pages <- pdf_text(pdf_path)
#     full_text <- paste(text_pages, collapse = "\n")
# 
#     # Save to file
#     writeLines(full_text, txt_path, useBytes = TRUE)
# 
#     cat("  ✅ Saved as:", basename(txt_path), "\n\n")
# 
#   }, error = function(e) {
#     cat("  ❌ Failed to process", filename, "–", e$message, "\n\n")
#   })
# }
# 
# cat("🎉 All done! Extracted texts saved in:\n", output_folder, "\n")
