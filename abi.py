
#!/usr/bin/python
# Copyright (C) 2014 Brian J. Stucky
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Script by Brian J. Stucky was slightly modified by Martin KM Engqvist



import SequenceTrace

class TraceFileError(Exception):
    pass

class UnknownFileTypeError(TraceFileError):
    def __str__(self):
        return 'The file format was not recognized.  Please convert the file to a supported format and try again.'


class ABIError(TraceFileError):
    pass

class ABIVersionError(ABIError):
    def __init__(self, ver_major, ver_minor):
        self.ver_major = ver_major
        self.ver_minor = ver_minor

    def __str__(self):
        return 'This file uses version ' + str(self.ver_major) + '.' + str(self.ver_minor) + ' of the ABI format.  This software only supports version 1.x of the format.'

class ABIIndexError(ABIError):
    def __init__(self, indexnum, indextotal):
        self.indexnum = indexnum
        self.indextotal = indextotal

    def __str__(self):
        return 'Error reading ABI file index entry ' + str(self.indexnum) + ' of ' + str(self.indextotal) + ' expected entries.  The file might be damaged.'

class ABIDataError(ABIError):
    def __init__(self, expectedlen, actuallen):
        self.expectedlen = expectedlen
        self.actuallen = actuallen

    def __str__(self):
        return 'Error reading ABI file data.  Expected ' + str(self.expectedlen) + ' bytes but only got ' + str(self.actuallen) + ' bytes.  The file appears to be damaged.'


class ABISequenceTrace(SequenceTrace):
    def loadFile(self, filename):
        self.fname = filename

        try:
            self.tf = open(filename, 'rb')
        except IOError:
            raise

        self.abiindex = list()

        # read the ABI magic number
        abinum = self.tf.read(4)
        #print abinum
        if abinum != 'ABIF':
            raise ABIError('The ABI file header is invalid.  The file appears to be damaged.')

        # check the major version number
        try:
            version = unpack('>H', self.tf.read(2))[0]
        except struct.error:
            raise ABIError('The ABI file header is invalid.  The file appears to be damaged.')
        #print version
        if (version / 100) != 1:
            raise ABIVersionError(version / 100, version % 100)
        
        # skip the next 10 bytes
        self.tf.read(10)
        
        # get the file index information
        try:
            index_entry_len = unpack('>h', self.tf.read(2))[0]
            self.num_index_entries = unpack('>i', self.tf.read(4))[0]
            total_index_size = unpack('>i', self.tf.read(4))[0]
            self.index_offset = unpack ('>i', self.tf.read(4))[0]
        except struct.error:
            raise ABIError('The ABI file header is invalid.  The file appears to be damaged.')
        
        #print index_entry_len, self.num_index_entries, total_index_size, self.index_offset

        self.readABIIndex()
        
        self.readBaseCalls()
        self.readConfScores()
        self.readTraceData()
        self.readBaseLocations()
        self.readComments()
        
    def readABIIndex(self):
        # read the ABI index block
        self.tf.seek(self.index_offset, 0)

        for cnt in range(self.num_index_entries):
            try:
                self.abiindex.append(dict(did=0, idv=0, dformat=0, fsize=0, dcnt=0, dlen=0, offset=0))
                self.abiindex[cnt]['did'] = self.tf.read(4)
                self.abiindex[cnt]['idv'] = unpack('>I', self.tf.read(4))[0]
                self.abiindex[cnt]['dformat'] = unpack('>H', self.tf.read(2))[0]
                self.abiindex[cnt]['fsize'] = unpack('>H', self.tf.read(2))[0]
                self.abiindex[cnt]['dcnt'] = unpack('>I', self.tf.read(4))[0]
                self.abiindex[cnt]['dlen'] = unpack('>I', self.tf.read(4))[0]
                self.abiindex[cnt]['offset'] = unpack('>I', self.tf.read(4))[0]
                # skip 4 bytes (the unused "data handle" field)
                self.tf.read(4)
            except struct.error:
                raise ABIIndexError(cnt, self.num_index_entries)
        
        #self.printABIIndex('CMNT')

    def printABIIndex(self, data_id):
        for entry in self.abiindex:
            if entry['did'] == data_id:
                print 'entry ID:', entry['did']
                print 'idv:', entry['idv']
                print 'data format:', entry['dformat']
                print 'format size:', entry['fsize']
                print 'data count:', entry['dcnt']
                print 'total data length:', entry['dlen']
                print 'data offset:', entry['offset']

    def getIndexEntry(self, data_id, number):
        for row in self.abiindex:
            if (row['did'] == data_id) and (row['idv'] == number):
                return row
    
        return None
    
    def getIndexEntriesById(self, data_id):
        entries = list()
    
        for row in self.abiindex:
            if row['did'] == data_id:
                entries.append(row)
    
        return entries
    
    # Attempts to get a bunch of information about the sequencing run from the ABI file.  As much as possible,
    # the keys used for individual comment values correspond with the keys used for the same values by the
    # Staden software package.  However, this method also retrieves some comments that are not read by the
    # Staden package.  To avoid confusion, these additional comment values are not given 4-letter keys.
    def readComments(self):
        # get the sample name
        entry = self.getIndexEntry('SMPL', 1)
        if entry:
            self.comments['NAME'] = self.readString(entry)

        # get the run name
        entry = self.getIndexEntry('RunN', 1)
        if entry:
            self.comments['Run name'] = self.readString(entry)

        # get the lane number
        entry = self.getIndexEntry('LANE', 1)
        if entry:
            self.comments['LANE'] = str(self.read2ByteInts(entry)[0])

        # get the signal strengths for each dye
        entry = self.getIndexEntry('S/N%', 1)
        if entry:
            stvals = self.read2ByteInts(entry)

            # use the "filter wheel order" to determine the base/value pairings
            order = self.getBaseDataOrder()
            sigst = {}
            for cnt in range(0, len(order)):
                sigst[order[cnt]] = stvals[cnt]

            self.comments['SIGN'] = 'A={0},C={1},G={2},T={3}'.format(sigst['A'], sigst['C'], sigst['G'], sigst['T'])

        # get the average peak spacing
        entry = self.getIndexEntry('SPAC', 1)
        if entry:
            spacing = self.read4ByteFloats(entry)[0]
            # if spacing is invalid, estimate it ourselves (the Staden code [seqIOABI.c] indicates this is a possibility)
            if spacing < 0:
                spacing = float(self.basepos[-1] - self.basepos[0]) / (len(self.basepos) - 1)
            self.comments['SPAC'] = '{0:.2f}'.format(spacing)

        # get the run dates and times
        d_entries = self.getIndexEntriesById('RUND')
        t_entries = self.getIndexEntriesById('RUNT')
        if (len(d_entries) > 1) and (len(t_entries) > 1):
            sdate = self.readDateTime(self.getIndexEntry('RUND', 1), self.getIndexEntry('RUNT', 1))
            edate = self.readDateTime(self.getIndexEntry('RUND', 2), self.getIndexEntry('RUNT', 2))
            #print sdate, edate
            self.comments['RUND'] = sdate.strftime('%Y%m%d.%H%M%S') + ' - ' + edate.strftime('%Y%m%d.%H%M%S')
            self.comments['DATE'] = sdate.strftime('%a %d %b %H:%M:%S %Y') + ' to ' + edate.strftime('%a %d %b %H:%M:%S %Y')

        # get the data collection dates and times
        if (len(d_entries) == 4) and (len(t_entries) == 4):
            sdate = self.readDateTime(self.getIndexEntry('RUND', 3), self.getIndexEntry('RUNT', 3))
            edate = self.readDateTime(self.getIndexEntry('RUND', 4), self.getIndexEntry('RUNT', 4))
            #print sdate, edate
            self.comments['Data coll. dates/times'] = sdate.strftime('%a %d %b %H:%M:%S %Y') + ' to ' + edate.strftime('%a %d %b %H:%M:%S %Y')

        # get the dye set/primer (mobility) file
        entry = self.getIndexEntry('PDMF', 1)
        if entry:
            self.comments['DYEP'] = self.readString(entry)

        # get the sequencing machine name and serial number
        entry = self.getIndexEntry('MCHN', 1)
        if entry:
            self.comments['MACH'] = self.readString(entry)

        # get the sequencing machine model
        entry = self.getIndexEntry('MODL', 1)
        if entry:
            self.comments['MODL'] = self.readString(entry)

        # get the basecaller name
        entry = self.getIndexEntry('SPAC', 2)
        if entry:
            self.comments['BCAL'] = self.readString(entry)

        # get the data collection software version
        entry = self.getIndexEntry('SVER', 1)
        if entry:
            self.comments['VER1'] = self.readString(entry)

        # get the basecaller version
        entry = self.getIndexEntry('SVER', 2)
        if entry:
            self.comments['VER2'] = self.readString(entry)

        # get the plate size
        entry = self.getIndexEntry('PSZE', 1)
        if entry:
            self.comments['Plate size'] = str(self.read4ByteInts(entry)[0])

        # get the gel name
        # This is included here because it is read by the Staden package, but it does not appear to be
        # included in the modern ABIF documentation.
        entry = self.getIndexEntry('GELN', 1)
        if entry:
            self.comments['GELN'] = self.readString(entry)

        # get the instrument (matrix) file
        # This is included here because it is read by the Staden package, but it does not appear to be
        # included in the modern ABIF documentation.
        entry = self.getIndexEntry('MTXF', 1)
        if entry:
            self.comments['MTXF'] = self.readString(entry)

        # 'APrX' points to a long XML string with detailed information about the analysis protocol used
        #entry = self.getIndexEntry('APrX', 1)
        #if entry:
        #    self.readString(entry)

    def readDateTime(self, dateindexrow, timeindexrow):
        # date format:
        #   bits 31-16: year
        #   bits 15-8: month
        #   bits 7-0: day of month
        # time format:
        #   bits 31-24: hour
        #   bits 23-16: minutes
        #   bits 15-8: seconds
        datenum = self.read4ByteInts(dateindexrow)[0]
        timenum = self.read4ByteInts(timeindexrow)[0]
        dateobj = datetime(year=(datenum >> 16), month=((datenum >> 8) & 0xff), day=(datenum & 0xff),
                hour=(timenum >> 24), minute=((timenum >> 16) & 0xff), second=((timenum >> 8) & 0xff))

        return dateobj

    def readString(self, indexrow):
        if indexrow['fsize'] != 1:
            raise ABIError('Index entry contains an invalid format size for string data.')
        if indexrow['dformat'] not in (2, 18, 19):
            raise ABIError('Index entry contains an invalid data type for character data.')
    
        if indexrow['dlen'] <= 4:
            # The actual data are stored in the offset field of the index entry.  Because the offset
            # was read as an unsigned, big-endian integer, the bytes should be in the correct order for
            # the following bit shift operations.
            lst = list()
            for cnt in range(0, indexrow['dcnt']):
                val = (indexrow['offset'] >> ((3 - cnt) * 8)) & 0xff
                lst.append(chr(val))
    
            strval = ''.join(lst)
        else:
            # get the data from the file
            self.tf.seek(indexrow['offset'], 0)
            strval = self.tf.read(indexrow['dcnt'])
    
        if indexrow['dlen'] != len(strval):
            raise ABIDataError(indexrow['dlen'], len(strval))

        # If this is a Pascal-style string (format 18), then remove the first character (which specifies
        # the string length).  If this is a C-style string (format 19), then remove the trailing
        # null character.
        if indexrow['dformat'] == 18:
            strval = strval[1:]
        elif indexrow['dformat'] == 19:
            strval = strval[:-1]

        return strval
    
    def read1ByteInts(self, indexrow):
        if indexrow['fsize'] != 1:
            raise ABIError('Index entry contains an invalid format size for 1-byte integers.')
    
        # see if the data format is signed or unsigned
        if indexrow['dformat'] == 1:
            formatstr = 'B'
        elif indexrow['dformat'] == 2:
            formatstr = 'b'
        else:
            raise ABIError('Index entry contains an invalid data type ID for 1-byte integers.')

        lst = list()
    
        if indexrow['dlen'] <= 4:
            # The actual data are stored in the offset field of the index entry.  Because the offset
            # was read as an unsigned, big-endian integer, the bytes should be in the correct order for
            # the following bit shift operations.
            # First, repack the integer to deal with the possibility of signed integers (shift operations
            # would only return positive values).
            data = pack('>I', indexrow['offset'])
            for cnt in range(0, indexrow['dcnt']):
                val = unpack(formatstr, data[cnt:cnt+1])[0]
                lst.append(val)
        else:
            # get the data from the file
            self.tf.seek(indexrow['offset'], 0)
            for cnt in range(0, indexrow['dcnt']):
                lst.append(unpack(formatstr, self.tf.read(1))[0])
    
        if indexrow['dlen'] != len(lst):
            raise ABIDataError(indexrow['dlen'], len(lst))

        return lst
    
    def read2ByteInts(self, indexrow):
        if indexrow['fsize'] != 2:
            raise ABIError('Index entry contains an invalid format size for 2-byte integers.')

        # see if the data format is signed or unsigned
        if indexrow['dformat'] == 3:
            formatstr = '>H'
        elif indexrow['dformat'] == 4:
            formatstr = '>h'
        else:
            raise ABIError('Index entry contains an invalid data type ID for 2-byte integers.')
    
        lst = list()
    
        if indexrow['dlen'] <= 4:
            # The actual data are stored in the offset field of the index entry.  Because the offset
            # was read as an unsigned, big-endian integer, the bytes should be in the correct order for
            # the following operations.
            # First, repack the integer to deal with the possibility of signed integers (shift operations
            # would only return positive values).
            data = pack('>I', indexrow['offset'])
            for cnt in range(0, indexrow['dcnt']):
                val = unpack(formatstr, data[cnt*2:cnt*2+2])[0]
                lst.append(val)
        else:
            # get the data from the file
            self.tf.seek(indexrow['offset'], 0)
            for cnt in range(0, indexrow['dcnt']):
                lst.append(unpack(formatstr, self.tf.read(2))[0])
    
        if indexrow['dlen'] != (len(lst) * 2):
            raise ABIDataError(indexrow['dlen'], (len(lst) * 2))
    
        return lst
    
    def read4ByteInts(self, indexrow):
        if indexrow['fsize'] != 4:
            raise ABIError('Index entry contains an invalid format size for 4-byte integers.')
        if indexrow['dformat'] not in (5, 10, 11):
            raise ABIError('Index entry contains an invalid data type ID for 4-byte integers.')
    
        lst = list()
    
        if indexrow['dlen'] == 4:
            # The actual data are stored in the offset field of the index entry.  In the case of 4-byte
            # ints, the offset value is the data value.  It must be repacked, though, to reinterpret it
            # as a signed integer.
            data = pack('>I', indexrow['offset'])
            val = unpack('>i', data)[0]
            lst.append(val)
        else:
            # get the data from the file
            self.tf.seek(indexrow['offset'], 0)
            for cnt in range(0, indexrow['dcnt']):
                lst.append(unpack('>i', self.tf.read(4))[0])
    
        if indexrow['dlen'] != (len(lst) * 4):
            raise ABIDataError(indexrow['dlen'], (len(lst) * 4))
    
        return lst

    def read4ByteFloats(self, indexrow):
        if indexrow['fsize'] != 4:
            raise ABIError('Index entry contains an invalid format size for 4-byte floating point numbers.')
        if indexrow['dformat'] != 7:
            raise ABIError('Index entry contains an invalid data type ID for 4-byte floating point numbers.')
    
        lst = list()
    
        if indexrow['dlen'] <= 4:
            # The actual data are stored in the offset field of the index entry.
            data = pack('>I', indexrow['offset'])
            lst.append(unpack('>f', data)[0])
        else:
            # get the data from the file
            self.tf.seek(indexrow['offset'], 0)
            for cnt in range(0, indexrow['dcnt']):
                lst.append(unpack('>f', self.tf.read(4))[0])
    
        if indexrow['dlen'] != (len(lst) * 4):
            raise ABIDataError(indexrow['dlen'], (len(lst) * 4))
    
        return lst

    # According to the ABIF documentation, ABIF files (after base calling) should contain two base
    # call entries ("PBAS"): one containing "sequence characters edited by user" (entry number 1),
    # and one containing "sequence characters as called by Basecaller" (entry number 2).  These
    # two entries will, in most cases, contain identical sequence data.  This method follows the same
    # convention used by the Staden package (see seqIOABI.c), which is to only look at entry 1 (the
    # user-edited sequence) and ignore entry 2.
    def readBaseCalls(self):
        row = self.getIndexEntry('PBAS', 1)
        if row is None:
            raise ABIError('No base call data were found in the ABI file.  The file might be damaged.')

        # read the base calls from the file
        self.basecalls = self.readString(row).upper()
    
    # There is an inconsistency in the ABIF file format documentation regarding the data format of the
    # confidence scores.  The data format ID (as actually found in a .ab1 file) is 2, indicating the
    # values are signed 1-byte integers.  The ABIF documentation, however, sugggests the values can range
    # from 0-255 (i.e., an unsigned 1-byte integer).  In practice, the actual values do not appear to
    # exceed 61, making the distinction between signed/unsigned irrelevant.  For now, the data format ID
    # is taken as the correct indication of the underlying data format.
    #
    # According to the ABIF documentation, ABIF files (after base calling) should contain two quality
    # value (QV) entries ("PCON"): one containing QVs "as edited by user" (entry number 1), and one
    # containing QVs "as called by Basecaller" (entry number 2).  These two entries will, in most cases,
    # contain identical values.  This method follows the same convention used by the Staden package
    # (see seqIOABI.c), which is to only look at entry 1 (the user-edited QVs) and ignore entry 2.
    def readConfScores(self):
        row = self.getIndexEntry('PCON', 1)
        if row is None:
            raise ABIError('No confidence score data were found in the ABI file.  SeqTrace requires confidence scores for all base calls.')
    
        # read the base call confidence scores from the file
        self.bcconf = self.read1ByteInts(row)
    
        return True
    
    # According to the ABIF documentation, ABIF files (after base calling) should contain two peak
    # location (PL) entries ("PLOC"): one containing PLs "edited by user" (entry number 1), and one
    # containing PLs "as called by Basecaller" (entry number 2).  These two entries will, in most cases,
    # contain identical information.  This method follows the same convention used by the Staden package
    # (see seqIOABI.c), which is to only look at entry 1 (the user-edited PLs) and ignore entry 2.
    def readBaseLocations(self):
        row = self.getIndexEntry('PLOC', 1)
        if row is None:
            raise ABIError('No base location data were found in the ABI file.  The file might be damaged.')
    
        # read the base call locations from the file
        self.basepos = self.read2ByteInts(row)

        return True

    def getBaseDataOrder(self):
        # retrieve the "filter wheel order" row from the file index
        rows = self.getIndexEntriesById('FWO_')
    
        if len(rows) > 1:
            raise ABIError('Found multiple filter wheel order index entries in ABI file.')
        if rows[0]['dlen'] != 4:
            raise ABIError('Incorrect data length for filter wheel order index entry.')
    
        # the data length is only 4 bytes, so the actual data is stored in the offset
        val = rows[0]['offset']
    
        base_order = list()
    
        base_order.append(chr((val >> 24) & 0xff))
        base_order.append(chr((val >> 16) & 0xff))
        base_order.append(chr((val >> 8) & 0xff))
        base_order.append(chr(val & 0xff))
    
        return base_order
    
    def readTraceData(self):
        base_order = self.getBaseDataOrder()
        maxval = 0
        
        # This is the ID for the first 'DATA' index entry that points to the processed
        # trace data.  The man page for the Staden program convert_trace states
        # that IDs 9-12 contain the processed data; IDs 1-4 contain the raw data.
        # The ABIF documentation from ABI also suggests that IDs 1-8 will always contain
        # raw data, and 9-12 will contain the processed data.  Is this always correct?
        start_id = 9
    
        for cnt in range(0, 4):
            row = self.getIndexEntry('DATA', start_id + cnt)
            if row == None:
                raise ABIError('Could not find trace data index entries for all bases.  The file might be damaged.')
    
            # read the trace data from the file
            lst = self.read2ByteInts(row)
            tmpmax = max(lst)
            if tmpmax > maxval:
                maxval = tmpmax
            self.tracesamps[base_order[cnt]] = lst

        self.max_traceval = maxval
