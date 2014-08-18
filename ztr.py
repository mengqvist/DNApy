
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


class ZTRError(TraceFileError):
    pass

class ZTRVersionError(ZTRError):
    def __init__(self, ver_major, ver_minor):
        self.ver_major = ver_major
        self.ver_minor = ver_minor

    def __str__(self):
        return 'This file uses version ' + str(self.ver_major) + '.' + str(self.ver_minor) + ' of the ZTR format.  This software only supports version 1.2 of the format.'

class ZTRDataFormatError(ZTRError):
    def __init__(self, format_id):
        self.format_id = format_id

    def __str__(self):
        return 'The ZTR data format ID ' + str(self.format_id) + ' is invalid or not supported.'

class ZTRMissingDataError(ZTRError):
    def __init__(self, expectedlen, actuallen):
        self.expectedlen = expectedlen
        self.actuallen = actuallen

    def __str__(self):
        return 'Error reading ZTR data chunk.  Expected ' + str(self.expectedlen) + ' bytes but only got ' + str(self.actuallen) + ' bytes.  The file appears to be damaged.'


class ZTRSequenceTrace(SequenceTrace):
    def loadFile(self, filename):
        self.fname = filename

        try:
            tf = open(filename, 'rb')
        except IOError:
            raise

        # read the header
        self.magicnum = tf.read(8)
        if self.magicnum != '\256ZTR\r\n\032\n':
            raise ZTRError('The ZTR file header is invalid.  The file appears to be damaged.')

        try:
            self.ver_major = unpack('b', tf.read(1))[0]
            self.ver_minor = unpack('b', tf.read(1))[0]
            #print 'major version number:', ver_major
            #print 'minor version number:', ver_minor
        except struct.error:
            raise ZTRError('The ZTR file header is invalid.  The file appears to be damaged.')

        if (self.ver_major != 1) or (self.ver_minor != 2):
            raise ZTRVersionError(self.ver_major, self.ver_minor)

        # read and process the data chunks
        chunk = self.readChunk(tf)
        while chunk != False:
            #print 'chunk type:', chunk[0]
            #print 'compressed data length:', chunk[1]
            #print 'uncompressed data length:', len(chunk[2])
            if chunk[0] == 'SMP4':
                # trace sample data
                self.readTraceSamples(chunk[2])
            elif chunk[0] == 'BASE':
                # base calls
                self.basecalls = chunk[2][1:].upper()
            elif chunk[0] == 'BPOS':
                # positions of base calls relative to trace samples
                self.basepos = list()
                # skip 4 leading null bytes
                for cnt in range(4, len(chunk[2]), 4):
                    self.basepos.append(unpack('>I', chunk[2][cnt:cnt+4])[0])
            elif chunk[0] == 'CNF4':
                # confidence scores; this is required to come after a BASE chunk
                self.bcconf = list()
                for cnt in range(1, self.getNumBaseCalls()+1):
                    self.bcconf.append(unpack('b', chunk[2][cnt])[0])
            elif chunk[0] == 'TEXT':
                # get the comment key/value strings, ignoring leading/trailing null characters
                keyvals = chunk[2][1:-2].split('\0')
                for cnt in range(0, len(keyvals), 2):
                    #print keyvals[cnt] + ': ' + keyvals[cnt+1]
                    self.comments[keyvals[cnt]] = keyvals[cnt+1]

            chunk = self.readChunk(tf)

    def readTraceSamples(self, chunkdata):
        self.max_traceval = 0
        tracelen = (len(chunkdata) - 2) / 4

        offset = 2
        basenum = 0
        for base in ['A','C','G','T']:
            thisbase = list()
            start = basenum*tracelen + offset
            for cnt in range(0, tracelen, 2):
                val = unpack('>H', chunkdata[start+cnt:start+cnt+2])[0]
                thisbase.append(val)

            tmpmax = max(thisbase)
            if tmpmax > self.max_traceval:
                self.max_traceval = tmpmax
            self.tracesamps[base] = thisbase
            basenum += 1

    def zlibUncompress(self, cdata):
        # In examining the Staden package source code, it appears that the 4-byte data length integer is not
        # guaranteed to be in big-endian byte order.  Might this be a bug in the Staden code?  For this reason,
        # native byte order is used here (in fact, using big-endian order will cause this to fail on an
        # x86 machine).
        udatalen = unpack('I', cdata[:4])[0]
    
        udata = zlib.decompress(cdata[4:])

        #print 'expected uncompressed data length:', udatalen
        #print 'actual uncompressed data length:', len(udata)
        if udatalen != len(udata):
            raise ZTRError('Zlib decompression failed.  The expected data length did not match the actual data length.')
        
        return udata
    
    def RLEUncompress(self, cdata):
        # In examining the Staden package source code, it appears that the 4-byte data length integer is not
        # guaranteed to be in big-endian byte order.  Might this be a bug in the Staden code?  For this reason,
        # native byte order is used here (in fact, using big-endian order will cause this to fail on an
        # x86 machine).
        udatalen = unpack('I', cdata[:4])[0]
        guard = cdata[4]
        #print unpack_from('=bIb', data[:6])
        #print 'guard byte:', guard
    
        cnt = 5
        udata = list()
        while cnt < len(cdata):
            if cdata[cnt] == guard:
                runlen = unpack('B', cdata[cnt+1])[0]
                #print 'run length:', runlen
                if runlen == 0:
                    udata.append(guard)
                    cnt += 2
                else:
                    runchar = cdata[cnt+2]
                    udata.extend(runlen*list(runchar))
                    cnt += 3
            else:
                udata.append(cdata[cnt])
                cnt += 1
    
        #print 'expected uncompressed data length:', udatalen
        #print 'actual uncompressed data length:', len(udata)
        if udatalen != len(udata):
            raise ZTRError('RLE decompression failed.  The expected data length did not match the actual data length.')
    
        return ''.join(udata)
    
    def followDecode(self, cdata):
        # read the decode table
        table = unpack('256B', cdata[:256])
        #print table
    
        udata = list()
        prev = unpack('B', cdata[256])[0]
        udata.append(cdata[256])
        for cnt in xrange(257, len(cdata)):
            diff = unpack('b', cdata[cnt])[0]
            actual = table[prev] - diff

            # simulate 1-byte unsigned overflow/underflow, if needed
            if actual < 0:
                actual += 256
            elif actual > 255:
                #print actual,(actual-256), table[prev], diff
                actual -= 256

            prev = actual
            #print actual
            udata.append(pack('B', actual))
        
        return ''.join(udata)
    
    def decode16To8(self, cdata):
        cnt = 0
        udata = list()
        while cnt < len(cdata):
            val = unpack('b', cdata[cnt])[0]
            if (val > -128) and (val < 128):
                udata.append(pack('>h', val))
                cnt += 1
            elif val == -128:
                #print unpack('b', cdata[cnt+1])[0]
                #print unpack('b', cdata[cnt+2])[0]
                udata.extend(cdata[cnt+1:cnt+3])
                cnt += 3
            else:
                raise ZTRError('Invalid value encountered while attempting to read 16- to 8-bit encoded ZTR data.')
    
        return ''.join(udata)
    
    def decode32To8(self, cdata):
        cnt = 0
        udata = list()
        while cnt < len(cdata):
            val = unpack('b', cdata[cnt])[0]
            if (val > -128) and (val < 128):
                udata.append(pack('>i', val))
                cnt += 1
            elif val == -128:
                #print unpack('b', cdata[cnt+1])[0]
                #print unpack('b', cdata[cnt+2])[0]
                udata.extend(cdata[cnt+1:cnt+5])
                cnt += 5
            else:
                raise ZTRError('Invalid value encountered while attempting to read 32- to 8-bit encoded ZTR data.')
    
        return ''.join(udata)
    
    def decode8BitDelta(self, cdata):
        levels = unpack('b', cdata[0])[0]
        #print 'levels:', levels
    
        # first, unpack the 1-byte values
        udata = list()
        for cnt in xrange(1, len(cdata)):
            val = unpack('B', cdata[cnt])[0]
            udata.append(val)

        # now apply the reverse delta filtering
        for clev in range(levels):
            prev = 0
            for cnt in xrange(0, len(udata)):
                actual = udata[cnt] + prev
                if actual > 255:
                    # simulate 1-byte integer overflow
                    actual -= 256
                prev = actual
                udata[cnt] = actual
    
        # repack the data
        tmpdata = list()
        for val in udata:
            tmpdata.append(pack('B', val))
    
        return ''.join(tmpdata)
    
    def decode16BitDelta(self, cdata):
        levels = unpack('b', cdata[0])[0]
        #print 'levels:', levels
    
        # first, unpack the 2-byte values
        udata = list()
        for cnt in xrange(1, len(cdata), 2):
            val = unpack('>H', cdata[cnt:cnt+2])[0]
            udata.append(val)

        # now apply the reverse delta filtering
        for clev in range(levels):
            prev = 0
            for cnt in xrange(0, len(udata)):
                actual = udata[cnt] + prev
                if actual > 65535:
                    # simulate 2-byte integer overflow
                    actual -= 65536
                prev = actual
                udata[cnt] = actual
    
        # repack the data
        tmpdata = list()
        for val in udata:
            tmpdata.append(pack('>H', val))
    
        return ''.join(tmpdata)
    
    def decode32BitDelta(self, cdata):
        levels = unpack('b', cdata[0])[0]
        #print 'levels:', levels
    
        # first, unpack the 4-byte values (skipping the 2 padding bytes)
        udata = list()
        for cnt in xrange(3, len(cdata), 4):
            val = unpack('>I', cdata[cnt:cnt+4])[0]
            udata.append(val)

        # now apply the reverse delta filtering
        for clev in range(levels):
            prev = 0
            for cnt in xrange(0, len(udata)):
                actual = udata[cnt] + prev
                if actual > 4294967295:
                    # simulate 1-byte integer overflow
                    actual -= 4294967296
                prev = actual
                udata[cnt] = actual
    
        # repack the data
        tmpdata = list()
        for val in udata:
            tmpdata.append(pack('>I', val))
    
        return ''.join(tmpdata)
    
    
    def readChunk(self, fp):
        # get the chunk descriptor
        chtype = fp.read(4)
        #print 'chunk type:', chtype

        # check for EOF
        if len(chtype) == 0:
            return False
        elif len(chtype) != 4:
            raise ZTRError('The ZTR data chunk type could not be read.  The file appears to be damaged.')
    
        try:
            mdlen = unpack('>I', fp.read(4))[0]
            #print 'metadata length:', mdlen
    
            # skip over the metadata
            fp.read(mdlen)
    
            # get the size of the data
            datalen = unpack('>I', fp.read(4))[0]
            #print 'data length:', datalen
        except struct.error:
            raise ZTRError('The ZTR data chunk header could not be read.  The file appears to be damaged.')

        # read the chunk data from the file
        data = fp.read(datalen)
        if datalen != len(data):
            raise ZTRMissingDataError(datalen, len(data))
    
        # iteratively process the chunk data until we get the "raw",
        # uncompressed data
        dataformat = unpack('b', data[0])[0]
        while dataformat != 0:
            #print 'data format:', dataformat
            if dataformat == 1:
                # run-length encoding
                data = self.RLEUncompress(data[1:])
            elif dataformat == 2:
                # zlib encoding
                data = self.zlibUncompress(data[1:])
            elif dataformat == 64:
                # 8-bit delta encoded
                data = self.decode8BitDelta(data[1:])
            elif dataformat == 65:
                # 16-bit delta encoded
                data = self.decode16BitDelta(data[1:])
            elif dataformat == 66:
                # 32-bit delta encoded
                data = self.decode32BitDelta(data[1:])
            elif dataformat == 70:
                # 16 to 8 bit conversion
                data = self.decode16To8(data[1:])
            elif dataformat == 71:
                # 32 to 8 bit conversion
                data = self.decode32To8(data[1:])
            elif dataformat == 72:
                # 'follow' encoding
                data = self.followDecode(data[1:])
            else:
                # invalid/unsupported data format
                raise ZTRDataFormatError(dataformat)

            dataformat = unpack('b', data[0])[0]
    
        return (chtype, datalen, data)
