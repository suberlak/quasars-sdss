import numpy
import os.path


def ReadColumns(FileName, ncols=None, delimiter='', comment_character='#', verbose=True):
    """ Read a number of columns from a file into numpy arrays.
    This routine tries to create arrays of floats--in this case 'nan' values
    are preserved. Columns that can't be cast as floats are returned as
    string arrays, with the caveat that all 'nan' datums in these string arrays
    are turned into Nones.
    
    ncols is the maximum number of columns to read
    delimiter, by default, is any whitespace
    if verbose is true, print out the number of lines read
    Raises IOError (for unreadable files or files with no readable lines)
    If you try to unpack too many columns, *you* will raise a ValueError
    """
    try:
        file = open(FileName)
    except IOError, e:
        if verbose:
            print "I can't open", FileName, "\nAborting..."
        raise IOError, "Can't open" + FileName + ": " + e
    nrows = 0
    rows = []
    for line in file:
        # if the line has no characters,
        # or the first non-whitespace character is a comment:
        if len(line.strip())==0 or line.strip()[0] == comment_character:
            continue
        line = line[:line.find(comment_character)] # return just the part before any '#'
        # split the line into 'words'
        if delimiter:
            line = line.split(delimiter)
        else:
            line = line.split()   # default is to split on any whitespace
        # if we only want some of these columns...
        if ncols:   
            line = line[:ncols]
        line = [SafeFloat(x) for x in line]  # carefully make a list of floats
        rows.append( line )         # rows is a list of lists
        nrows += 1
    file.close()
    if verbose:
        print "Read", nrows, "lines from", FileName
    if len(rows) == 0:
        # oops, no readable lines
        raise IOError, "No readable lines in file"
    cols = zip( *rows ) # return a list of each column
                        # ( the * unpacks the 1st level of the list )
    outlist = []
    for col in cols:
        arr = numpy.asarray( col )
        type = arr.dtype
        if str(type)[0:2] == '|S':
            # it's a string array!
            outlist.append( arr )
        else:
            outlist.append( numpy.asarray(arr, numpy.float32) ) # it's (hopefully) just numbers
    return outlist # a list of numpy arrays, one for each column


def WriteColumns(FileName, cols, header='', clobber=False, NaNstring=None, verbose=True):
    """ Write a file containing a column for each array passed.

    Takes a filename, the arrays (which should be passed as a tuple),
    an optional header, and options enabling overwriting files, dealing
    with null values, and verbosity (defaults: protect files, use
    literal 'nan's, and tell you stuff ).
    If you don't want literal 'nan's to appear in your output file,
    set NaNstring to the *string* you wish to see instead:

      NaNstring='-99'   <- note that '-99' is a string, not an integer

    The header text is inserted verbatim (plus a final carriage return),
    so be sure to include comment characters if necessary.
        
    Raises IOError with a hopefully helpful error string

    BUGS: unexpected results if you give it a single string array outside
    of a tuple (try it!)
    """
    # should we check if we're going to overwrite something?
    if clobber == False and os.path.exists(FileName): 
        raise IOError, FileName + " already exists!"
    try:
        rows = len(cols[0])
    except TypeError:
        # there was only one column,
        cols = (cols,) # so turn cols into the expected tuple
        rows = len(cols[0])
    try:
        File = open(FileName, "w")
        if header != '':
            File.write( header+"\n" )
        for i in range(rows):
            for array in cols:
                # if NaNstring is defined, replace NaNs
                if NaNstring:
                    datum = array[i]
                    if numpy.isnan(datum)==True: 
                        File.write( "\t" + NaNstring )
                    else:
                        File.write( "\t" + str(array[i]) )
                else: # or don't
                    File.write( "\t" + str(array[i]) )
            File.write( "\n" )
    except IOError, e:
        File.close()
        if verbose:
            print "Error writing to", FileName
        raise IOError, "Error writing to " + FileName + ": " + e
    if verbose:
        print FileName, ":", rows, "rows written."
    File.close()


def SafeFloat(value):
    """ Try to convert values to floats, but recognize 'nan' and giving those
    back as None. Strings that can't be converted to floats are simply returned
    Used by ReadColumns()
    """
    try:
        return float(value)
    except:
        if value == 'nan':
            return None
        return value

