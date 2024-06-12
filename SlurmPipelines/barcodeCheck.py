import os

def extract_barcode_from_header(header, barcode_length):
    """Extracts barcode sequence from FASTQ read header."""
    return header.split(":")[-1][:barcode_length]

def check_barcode(barcode, expected_barcodes, allowed_mismatches=0):
    """Checks if the given barcode matches any of the expected barcodes."""
    for expected_barcode in expected_barcodes:
        mismatch_count = sum(1 for a, b in zip(barcode, expected_barcode) if a != b)
        if mismatch_count <= allowed_mismatches:
            return True
    return False

def process_fastq_pair(fastq_file_r1, fastq_file_r2, barcode_length, expected_barcodes, output_folder):
    """Process a pair of FASTQ files (R1 and R2)."""
    output_file_r1 = os.path.join(output_folder, os.path.splitext(os.path.basename(fastq_file_r1))[0] + "_barcode_check_R1.txt")
    output_file_r2 = os.path.join(output_folder, os.path.splitext(os.path.basename(fastq_file_r2))[0] + "_barcode_check_R2.txt")
    
    with open(fastq_file_r1, 'r') as f1, open(fastq_file_r2, 'r') as f2, \
         open(output_file_r1, 'w') as out_r1, open(output_file_r2, 'w') as out_r2:
        
        out_r1.write("Barcode Check Results for {}\n".format(os.path.basename(fastq_file_r1)))
        out_r2.write("Barcode Check Results for {}\n".format(os.path.basename(fastq_file_r2)))
        
        line_number_r1 = 0
        line_number_r2 = 0
        current_sequence_r1 = ''
        current_sequence_r2 = ''
        current_header_r1 = ''
        current_header_r2 = ''
        current_quality_r1 = ''
        current_quality_r2 = ''
        
        for line_r1, line_r2 in zip(f1, f2):
            line_number_r1 += 1
            line_number_r2 += 1
            
            if line_number_r1 % 4 == 1:
                current_header_r1 = line_r1.strip()
            elif line_number_r1 % 4 == 2:
                current_sequence_r1 = line_r1.strip()
            elif line_number_r1 % 4 == 0:
                current_quality_r1 = line_r1.strip()
                
                # Extract barcode from R1 header
                barcode_r1 = extract_barcode_from_header(current_header_r1, barcode_length)
                
                # Extract barcode from R2 header
                barcode_r2 = extract_barcode_from_header(current_header_r2, barcode_length)
                
                # Check if barcodes match expected barcodes
                if check_barcode(barcode_r1, expected_barcodes) and check_barcode(barcode_r2, expected_barcodes):
                    out_r1.write("Barcode: {}, Sequence: {}\n".format(barcode_r1, current_sequence_r1))
                    out_r2.write("Barcode: {}, Sequence: {}\n".format(barcode_r2, current_sequence_r2))
                else:
                    out_r1.write("Mismatched Barcode: {}, Sequence: {}\n".format(barcode_r1, current_sequence_r1))
                    out_r2.write("Mismatched Barcode: {}, Sequence: {}\n".format(barcode_r2, current_sequence_r2))

# Process all FASTQ pairs in a folder
def process_fastq_folder(folder_path, barcode_length, expected_barcodes, output_folder):
    """Process all FASTQ files in a folder."""
    for file_name in os.listdir(folder_path):
        if file_name.endswith("_R1.fastq") or file_name.endswith("_R1.fq"):
            fastq_file_r1 = os.path.join(folder_path, file_name)
            fastq_file_r2 = os.path.join(folder_path, file_name.replace("_R1", "_R2"))
            if os.path.exists(fastq_file_r2):
                process_fastq_pair(fastq_file_r1, fastq_file_r2, barcode_length, expected_barcodes, output_folder)

# Example usage
if __name__ == "__main__":
    folder_path = "./"
    output_folder = "./"
    barcode_length = 6
    expected_barcodes = ["ATCGAA", "CGATGC", "TACGTA"]  # Replace with your expected barcode sequences
    process_fastq_folder(folder_path, barcode_length, expected_barcodes, output_folder)
