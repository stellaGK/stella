
path = 'default_stella_input.in'

# Open text file
f = open(path, "r+")

# Read the text file
text = f.read()
lines = text.split('\n')

# Find the knobs
knobs = []
for line in lines:
    if len(line)>0:
        if line[0]=='&':
            knobs.append(line)
            
# Write the knobs in lower case
for knob in knobs:
    text_before = text.split(knob)[0]
    text_after = f'{knob}'.join(text.split(knob)[1:])
    text_after = '\n /'.join(text_after.split('\n /')[1:])
    text_knob = text.split(knob)[1].split('\n /')[0] 
    text_knob = text_knob.lower()
    text_knob = text_knob.replace('\n', '\n ')
    text = text_before + knob.lower() + text_knob + '\n /' + text_after

# Make the text pretty
text = text.replace('=', ' = ')
text = text.replace(' = f,', ' = .false.')
text = text.replace(' = t,', ' = .true.')

# Remove useless spaces within quotes
lines = text.split('\n')
for i, line in enumerate(lines):
    if '"' in line:
        line = line.split('"')[0] + '"' + line.split('"')[1].replace(' ','') + '"' + line.split('"')[2].replace(',','')
        lines[i] = line
text = '\n'.join(lines)
text = text.replace('\n   ', '\n  ')
text = text.replace('\n   ', '\n  ') 

# Clean up floats 
text = text.replace(' =  ', ' = ')
text = text.replace(' =  ', ' = ')
text = text.replace(' =  ', ' = ')
lines = text.split('\n')
for i, line in enumerate(lines):
    if ',' in line and ' = ' in line and '<' not in line: 
        before = line.split(' = ')[0] 
        number = line.split(' = ')[1].split(',')[0]
        number = float(number)
        if round(number,6)==number: number = round(number,5)
        line = before + ' = ' + str(number)
        lines[i] = line
text = '\n'.join(lines)

# Write to text file
f.seek(0)
f.write(text)
f.truncate() 
print(f'\nCleaned up the "{path}" file.\n')

# Close file
f.close()


