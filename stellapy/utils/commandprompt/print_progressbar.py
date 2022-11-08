
# Print iterations progress
def print_progressbar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    if len(suffix)<50: suffix = suffix+(50-len(suffix))*" "
    if iteration==0: print(f'\n\r{prefix} |{bar}| {percent}% {suffix} \n', end = printEnd) 
    if iteration!=0: print(f'\033[A\r{prefix} |{bar}| {percent}% {suffix} \n', end = printEnd) 
    # Print New Line on Complete
    if iteration == total: 
        print()