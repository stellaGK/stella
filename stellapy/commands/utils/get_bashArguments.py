
import sys, getopt
def get_bashArguments(options, long_options, print_help, function, default_arguments):    
    ''' Read the options and arguments from the terminal.
    
    opts_args, comments = getopt.getopt(args, options, [long_options])
    '''
    try:
        opts_args = getopt.getopt(sys.argv[1:],options,long_options)[0]
        return opts_args
    except getopt.GetoptError:
        print('Internal error')
        if print_help:
            print_help(function, default_arguments)
        sys.exit(2)