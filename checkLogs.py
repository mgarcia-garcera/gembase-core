def parse_logfile(filepath):
    input_file = None
    success = False
    error_code = ""
    error_summary = ""

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        for line in lines:
            if line.startswith("my input is"):
                input_file = line.strip().split("is", 1)[1].strip()
        
        # Check success in the last line
        if lines and "The pipeline was successful" in lines[-1]:
            success = True
        else:
            # Check for known errors
            for line in lines:
                if "Log file not found" in line:
                    error_code = "E001"
                    error_summary = "Missing internal log file"
                    break
                if "Pipeline did not complete successfully" in line:
                    error_code = "E002"
                    error_summary = "Pipeline failure or missing log"
                    break
            if not error_code:
                error_code = "E999"
                error_summary = "Unknown error"
    except Exception as e:
        input_file = filepath
        error_code = "E900"
        error_summary = f"Error reading file: {str(e)}"

    if success:
        return f"{input_file}\tOK\tPipeline completed successfully"
    else:
        return f"{input_file}\t{error_code}\t{error_summary}"


def main(listfile_path):
    try:
        with open(listfile_path, 'r') as f:
            log_files = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Failed to read listfile: {e}")
        return

    for log_file in log_files:
        result = parse_logfile(log_file)
        print(result)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python check_logs.py <listfile.txt>")
    else:
        main(sys.argv[1])