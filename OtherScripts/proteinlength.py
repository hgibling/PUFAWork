lengths = dict()
with open('fastainfo2.txt') as f:
    for line in f:
        if line.startswith('>'):
            #line.rstrip()
            lengthpos = line.rfind(':')
            aalength = line[lengthpos+1:]
            aalength = aalength.rstrip()
            descpos1 = line.find('|')
            newdesc = line[descpos1+1:]
            descpos2 = newdesc.find('|')
            desc = newdesc[:descpos2]
            lengths[desc] = aalength
   
maxID = max(lengths, key=lengths.get)

saved = open('finalproteins.txt', 'a')
saved.write(maxID)
saved.write('\n')