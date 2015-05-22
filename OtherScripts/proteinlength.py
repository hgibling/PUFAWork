with open('fastainfo.txt') as f:
    for line in f:
        if line.startswith('>'):
            lengthpos = line.rfind(':')
            print line[lengthpos+1:]
        
