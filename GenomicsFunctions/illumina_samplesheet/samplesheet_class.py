#!/usr/bin/env python

# Objects to read and write Illumina SampleSheets.
# SampleSheetReader(infilename)
# SampleSheetWriter(SampleSheetReader, outfilename, <column seperator (default ',')>, number of seperators)

import collections


class SampleSheetReader(object):
    """
    Read an Illumina SampleSheet and parse in to usable data structure.
    """

    def __init__(self, filename, sep=','):

        self.header_titles = ("Header", "Reads", "Settings", "Data")
        self.filename = filename

        self._sep = sep  # Define separator, usually comma.
        self.sep_count = self._val_count_cols()  # Check to make sure all rows have same number of columns.

        self.header = self._populate_sections(filename, self.header_titles[0], sep) 
        self.runid = self.header['Experiment Name'][0]
        self.date = self.header['Date'][0]
        self.assay = self.header['Assay'][0]
        
        self.reads = self._populate_sections(filename, self.header_titles[1], sep) 
        self.settings = self._populate_sections(filename, self.header_titles[2], sep) 

        # Special case of Data section.
        self.data = self._populate_sections(filename, self.header_titles[3], sep) 
        self._data_header = self._create_data_header()
        self.data = self._reformat_data_section(self.data, self._data_header)

    def __iter__(self):
        """
        Define iteration over SampleSheet field dictionaries.
        """

        field_list = (self.header, self.reads, self.settings, self.data)

        for field in field_list:
            yield field

    def _create_data_header(self):
        """
        Return a list that is the header of the Data section.
        """

        try:
            data_header = self.data['Sample_ID'][0]
        except:
            raise KeyError("Key Sample_ID is not contained within the SampleSheet Data section.")

        return data_header

    def _reformat_data_section(self, data_dict, data_header):
        """
        Break the ordered dictionary up in to a new data structure that facilitates calling individual values from this
        section.
        data_header = [str1, str2, ...]
        data_dict = {'SampleID': [[data]]}
        """
        
        sample_dict = collections.OrderedDict()
        for key in self.data:
            if key != self._data_header[0]:
                sample_dict[key] = collections.OrderedDict({self._data_header[0]: self.data[key][0][0]})
                for i in range(1, len(data_header)):
                    sample_dict[key][data_header[i]] = data_dict[key][0][i]

        return sample_dict

    def _val_count_cols(self):
        """
        Check to make sure each row contains the same number of columns (defined by sep).
        """

        first = True
 
        handle = open(self.filename, 'rU')
       
        with handle as samplesheet:
            for line in samplesheet:

                if first:
                    first_count = line.count(self._sep)
                    first = False
                else:
                    if line.count(self._sep) != first_count:
                        raise Exception("SampleSheet malformed, column count differs between rows.")
                    
        return first_count

    def _populate_sections(self, filename, header_title, sep):
        """
        Fill each section of the SampleSheet in to a dictionary.
        """
        
        local_dict = collections.OrderedDict()
        filling = False

        reader = open(filename, 'rU')

        with reader as samplesheet:
            for line in samplesheet:
                if (filling and
                    line.startswith('[')
                    and line.replace(',', '').rstrip(
                    '\n')[1:-1] != header_title.title()):
                    filling = False
                elif filling:
                    line = line.rstrip('\n').split(sep)
                    if line[0] not in local_dict:
                        local_dict[line[0]] = [line]
                    else:
                        local_dict[line[0]].append(line)
                elif not filling and \
                        line.startswith('[') and \
                        (line.replace(',', '')
                                 .rstrip('\n')[1:-1] == header_title.title()):
                    filling = True

        return local_dict


class SampleSheetWriter(object):
    """
    Write an Illumina SampleSheet.
    """

    def __init__(self, samplesheet, filename, sep=','):

        self.samplesheet = samplesheet
        self.filename = filename
        self.sep = sep
        self.write_sheet(filename, samplesheet, sep, samplesheet.sep_count)


    def write_sheet(self, filename, samplesheet, sep, sep_count):
        """
        Write the new SampleSheet from the samplesheet object.
        """

        handle_out = open(filename, 'w')

        i = 0
        no_data_header = True
        for field in samplesheet:
            handle_out.write(self._prepare_section_header(samplesheet.header_titles[i]))
            handle_out.write(sep*sep_count)
            handle_out.write('\n')
            if field != samplesheet.data:
                for entry in field:
                    print(entry)
                    print(field[entry][0])
                    if len(field[entry]) == 2:
                        for duplicate in field[entry]:
                            handle_out.write(sep.join(duplicate))
                            handle_out.write(self._fill_seps(field[entry][0],
                                                             sep_count))
                            handle_out.write('\n')
                    else:
                        handle_out.write(sep.join(field[entry][0]))
                        handle_out.write(self._fill_seps(field[entry][0],
                                                         sep_count))
                        handle_out.write('\n')
            else:
                for sample in field:

                    if no_data_header:
                        handle_out.write(sep.join([entry for entry in field[sample]]))
                        handle_out.write('\n')
                        no_data_header = False

                    if not no_data_header:
                        handle_out.write(sep.join([field[sample][entry] for entry in field[sample]]))
                        handle_out.write('\n')

            i += 1

        handle_out.close()

    def _fill_seps(self, entry, sep_count):
        """
        :return:
        """

        if (len(entry)-1) < sep_count:
            num_sep = sep_count - (len(entry)-1)
            return (num_sep*self.sep)

        return ""

    def _prepare_section_header(self, title):
        """
        Prepare a section title to be written.
        """
        new_title = '[' + title.title() + ']'
        return new_title
