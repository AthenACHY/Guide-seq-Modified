import svgwrite
import sys
import os

boxWidth = 10
box_size = 15
v_spacing = 3

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', "R": '#B3B3B3'}


def parseSitesFile(infile):
    offtargets = []
    with open(infile, 'r') as f:
        f.readline().rstrip()
        for line in f:
            line_items = line.split('\t')
            offtarget_sequence = line_items[21]
            offtarget_reads = line_items[11]
            ref_seq = line_items[32]
	    mismatches=line_items[22]
            if offtarget_sequence != '':
                offtargets.append({'seq': offtarget_sequence.strip(),
                                   'reads': int(offtarget_reads.strip()),
				   'mismatches': int(mismatches.strip())})
    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)
    return offtargets, ref_seq


def visualizeOfftargets(infile, outfile, title=None):

    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get offtargets array from file
    offtargets, ref_seq = parseSitesFile(infile)
    print ref_seq
    print offtargets
    # Initiate canvas
    dwg = svgwrite.Drawing(outfile + '.svg', profile='full')

    if title is not None:
        # Define top and left margins
        x_offset = 20
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 20
        y_offset = 20
# Draw ticks
    tick_locations = range(1, len(ref_seq) -5, 2)
    tick_locations.sort()

### odd no vs even no
    if len(ref_seq) %2==1:
        for x, y in zip(tick_locations[::-1], tick_locations):
            dwg.add(dwg.text(str(x), insert=(x_offset + (y - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))
    else:
        for x, y in zip(tick_locations[::-1], tick_locations):
            dwg.add(dwg.text(str(x), insert=(x_offset + y * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

	
# Draw reference sequence row
    for i, c in enumerate(ref_seq):
	print (i, c)
        y = y_offset
        x = x_offset + i * box_size
        try:
		dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        	dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))
        	dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(ref_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
	except:
		continue

    # Draw aligned sequence rows
    y_offset += 10  # leave some extra space after the reference row
    for j, seq in enumerate(offtargets):
        y = y_offset + j * box_size
        for i, c in enumerate(seq['seq']):
            x = x_offset + i * box_size
            if c == ref_seq[i] or ref_seq[i] == 'N':
                dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
            else:
                dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
        if seq['mismatches']==0:
		reads_text = dwg.text(str(seq['reads']) + "*", insert=(box_size * (len(ref_seq) + 1) + 20, y_offset + box_size * (j + 2) - 2), fill='black', style="font-size:15px; font-family:Courier")
        else:
                reads_text = dwg.text(str(seq['reads']), insert=(box_size * (len(ref_seq) + 1) + 20, y_offset + box_size * (j + 2) - 2), fill='black', style="font-size:15px; font-family:Courier")
        dwg.add(reads_text)

    dwg.save()


def main():
    if len(sys.argv) >= 3:
        if len(sys.argv) == 4:
            title = sys.argv[3]
        else:
            title = None
        visualizeOfftargets(sys.argv[1], sys.argv[2], title=title)
    else:
        print 'Usage: python visualization.py INFILE OUTFILE [TITLE]'
    from svglib.svglib import svg2rlg
    from reportlab.graphics import renderPDF
    import svgutils
    outfile=sys.argv[2]
    svg = svgutils.transform.fromfile(outfile + '.svg')
    originalSVG = svgutils.compose.SVG(outfile + '.svg')
    drawing = svg2rlg(outfile + ".svg")
    renderPDF.drawToFile(drawing, outfile + ".pdf", autoSize=0)
#    import cairosvg
#    cairosvg.svg2pdf(url=outfile +'.svg', write_to=outfile +'.pdf')

if __name__ == '__main__':
    main()

