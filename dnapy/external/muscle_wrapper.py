# Copyright 2001-2004 Brad Chapman. 
# Revisions copyright 2009-2013 by Peter Cock.
# Copyright 2009 by Cymon J. Cox.  All rights reserved. 
# This code is part of the Biopython distribution and governed by its 
# license.  Please see the LICENSE file that should have been included 
# as part of this package.


#                 Biopython License Agreement
#
#Permission to use, copy, modify, and distribute this software and its
#documentation with or without modifications and for any purpose and
#without fee is hereby granted, provided that any copyright notices
#appear in all copies and that both those copyright notices and this
#permission notice appear in supporting documentation, and that the
#names of the contributors or copyright holders not be used in
#advertising or publicity pertaining to distribution of the software
#without specific prior permission.
#
#THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
#WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
#CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
#OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
#OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
#OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
#OR PERFORMANCE OF THIS SOFTWARE.

#This file has been modified by Martin Engqvist for the DNApy project.

 
"""Command line wrapper for the multiple alignment program MUSCLE. 
""" 
 
 
 
__docformat__ = "restructuredtext en"  # Don't just use plain text in epydoc API pages! 
 
 
import os 
import platform 
import sys 
import subprocess 
import re 
from subprocess import CalledProcessError as _ProcessCalledError 
 



"""General mechanisms to access applications in Biopython. 
 
This module is not intended for direct use. It provides the basic objects which 
are subclassed by our command line wrappers, such as: 
 
 - Bio.Align.Applications 
 - Bio.Blast.Applications 
 - Bio.Emboss.Applications 
 - Bio.Sequencing.Applications 
 
These modules provide wrapper classes for command line tools to help you 
construct command line strings by setting the values of each parameter. 
The finished command line strings are then normally invoked via the built-in 
Python module subprocess. 
""" 
 
 
#Use this regular expression to test the property names are going to 
#be valid as Python properties or arguments 
_re_prop_name = re.compile(r"^[a-zA-Z][a-zA-Z0-9_]*$") 
assert _re_prop_name.match("t") 
assert _re_prop_name.match("test") 
assert _re_prop_name.match("_test") is None # we don't want private names 
assert _re_prop_name.match("-test") is None 
assert _re_prop_name.match("any-hyphen") is None 
assert _re_prop_name.match("underscore_ok") 
assert _re_prop_name.match("test_name") 
assert _re_prop_name.match("test2") 

#These are reserved names in Python itself, 
_reserved_names = ["and", "del", "from", "not", "while", "as", "elif", 
                   "global", "or", "with", "assert", "else", "if", "pass", 
                   "yield", "break", "except", "import", "print", "class", 
                   "exec", "in", "raise", "continue", "finally", "is", 
                   "return", "def", "for", "lambda", "try"] 

#These are reserved names due to the way the wrappers work 
_local_reserved_names = ["set_parameter"] 
 
 
class ApplicationError(_ProcessCalledError): 
    """Raised when an application returns a non-zero exit status. 
 
    The exit status will be stored in the returncode attribute, similarly 
    the command line string used in the cmd attribute, and (if captured) 
    stdout and stderr as strings. 
 
    This exception is a subclass of subprocess.CalledProcessError. 
 
    >>> err = ApplicationError(-11, "helloworld", "", "Some error text") 
    >>> err.returncode, err.cmd, err.stdout, err.stderr 
    (-11, 'helloworld', '', 'Some error text') 
    >>> print(err) 
    Command 'helloworld' returned non-zero exit status -11, 'Some error text' 
    """ 
    def __init__(self, returncode, cmd, stdout="", stderr=""): 
        self.returncode = returncode 
        self.cmd = cmd 
        self.stdout = stdout 
        self.stderr = stderr 
 
    def __str__(self): 
        #get first line of any stderr message 
        try: 
            msg = self.stderr.lstrip().split("\n", 1)[0].rstrip() 
        except: 
            msg = "" 
        if msg: 
            return "Command '%s' returned non-zero exit status %d, %r" % (self.cmd, self.returncode, msg) 
        else: 
            return "Command '%s' returned non-zero exit status %d" % (self.cmd, self.returncode) 
 
    def __repr__(self): 
        return "ApplicationError(%i, %s, %s, %s)" % (self.returncode, self.cmd, self.stdout, self.stderr) 
 
 
class AbstractCommandline(object): 
    """Generic interface for constructing command line strings. 
 
    This class shouldn't be called directly; it should be subclassed to 
    provide an implementation for a specific application. 
 
    For a usage example we'll show one of the EMBOSS wrappers.  You can set 
    options when creating the wrapper object using keyword arguments - or 
    later using their corresponding properties: 
 
    >>> from Bio.Emboss.Applications import WaterCommandline 
    >>> cline = WaterCommandline(gapopen=10, gapextend=0.5) 
    >>> cline 
    WaterCommandline(cmd='water', gapopen=10, gapextend=0.5) 
 
    You can instead manipulate the parameters via their properties, e.g. 
 
    >>> cline.gapopen 
    10 
    >>> cline.gapopen = 20 
    >>> cline 
    WaterCommandline(cmd='water', gapopen=20, gapextend=0.5) 
 
    You can clear a parameter you have already added by 'deleting' the 
    corresponding property: 
 
    >>> del cline.gapopen 
    >>> cline.gapopen 
    >>> cline 
    WaterCommandline(cmd='water', gapextend=0.5) 
 
    Once you have set the parameters you need, you can turn the object into 
    a string (e.g. to log the command): 
 
    >>> str(cline) 
    Traceback (most recent call last): 
    ... 
    ValueError: You must either set outfile (output filename), or enable filter or stdout (output to stdout). 
 
    In this case the wrapper knows certain arguments are required to construct 
    a valid command line for the tool.  For a complete example, 
 
    >>> from Bio.Emboss.Applications import WaterCommandline 
    >>> water_cmd = WaterCommandline(gapopen=10, gapextend=0.5) 
    >>> water_cmd.asequence = "asis:ACCCGGGCGCGGT" 
    >>> water_cmd.bsequence = "asis:ACCCGAGCGCGGT" 
    >>> water_cmd.outfile = "temp_water.txt" 
    >>> print(water_cmd) 
    water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5 
    >>> water_cmd 
    WaterCommandline(cmd='water', outfile='temp_water.txt', asequence='asis:ACCCGGGCGCGGT', bsequence='asis:ACCCGAGCGCGGT', gapopen=10, gapextend=0.5) 
 
    You would typically run the command line via a standard Python operating 
    system call using the subprocess module for full control. For the simple 
    case where you just want to run the command and get the output: 
 
    stdout, stderr = water_cmd() 
 
    Note that by default we assume the underlying tool is installed on the 
    system $PATH environment variable. This is normal under Linux/Unix, but 
    may need to be done manually under Windows. Alternatively, you can specify 
    the full path to the binary as the first argument (cmd): 
 
    >>> from Bio.Emboss.Applications import WaterCommandline 
    >>> water_cmd = WaterCommandline("C:\Program Files\EMBOSS\water.exe", 
    ...                              gapopen=10, gapextend=0.5, 
    ...                              asequence="asis:ACCCGGGCGCGGT", 
    ...                              bsequence="asis:ACCCGAGCGCGGT", 
    ...                              outfile="temp_water.txt") 
    >>> print(water_cmd) 
    "C:\Program Files\EMBOSS\water.exe" -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5 
 
    Notice that since the path name includes a space it has automatically 
    been quoted. 
 
    """ 
    #TODO - Replace the above example since EMBOSS doesn't work properly 
    #if installed into a folder with a space like "C:\Program Files\EMBOSS" 
 
    #Note the call example above is not a doctest as we can't handle EMBOSS 
    #(or any other tool) being missing in the unit tests. 
    def __init__(self, cmd, **kwargs): 
        """Create a new instance of a command line wrapper object.""" 
        # Init method - should be subclassed! 
        # 
        # The subclass methods should look like this: 
        # 
        # def __init__(self, cmd="muscle", **kwargs): 
        #     self.parameters = [...] 
        #     AbstractCommandline.__init__(self, cmd, **kwargs) 
        # 
        # i.e. There should have an optional argument "cmd" to set the location 
        # of the executable (with a sensible default which should work if the 
        # command is on the path on Unix), and keyword arguments.  It should 
        # then define a list of parameters, all objects derived from the base 
        # class _AbstractParameter. 
        # 
        # The keyword arguments should be any valid parameter name, and will 
        # be used to set the associated parameter. 
        self.program_name = cmd 
        try: 
            parameters = self.parameters 
        except AttributeError: 
            raise AttributeError("Subclass should have defined self.parameters") 
        #Create properties for each parameter at run time 
        aliases = set() 
        
        for p in parameters:

            if not p.names: 
                assert isinstance(p, _StaticArgument), p 
                continue 
            for name in p.names: 
                if name in aliases: 
                    raise ValueError("Parameter alias %s multiply defined" 
                                     % name) 
                aliases.add(name) 
            name = p.names[-1] 
            if _re_prop_name.match(name) is None: 
                raise ValueError("Final parameter name %s cannot be used as " 
                                 "an argument or property name in python" 
                                 % repr(name)) 
            if name in _reserved_names: 
                raise ValueError("Final parameter name %s cannot be used as " 
                                 "an argument or property name because it is " 
                                 "a reserved word in python" % repr(name)) 
            if name in _local_reserved_names: 
                raise ValueError("Final parameter name %s cannot be used as " 
                                 "an argument or property name due to the " 
                                 "way the AbstractCommandline class works" 
                                 % repr(name)) 
 
            #Beware of binding-versus-assignment confusion issues 
            def getter(name): 
                return lambda x: x._get_parameter(name) 
 
            def setter(name): 
                return lambda x, value: x.set_parameter(name, value) 
 
            def deleter(name): 
                return lambda x: x._clear_parameter(name) 
 
            doc = p.description 
            if isinstance(p, _Switch): 
                doc += "\n\nThis property controls the addition of the %s switch, treat this property as a boolean." % p.names[0] 
            else: 
                doc += "\n\nThis controls the addition of the %s parameter and its associated value.  Set this property to the argument value required." % p.names[0] 
            prop = property(getter(name), setter(name), deleter(name), doc) 
            setattr(self.__class__, name, prop)  # magic! 

        for key, value in list(kwargs.items()):
            self.set_parameter(key, value) 

 
    def _validate(self): 
        """Make sure the required parameters have been set (PRIVATE). 
 
        No return value - it either works or raises a ValueError. 
 
        This is a separate method (called from __str__) so that subclasses may 
        override it. 
        """ 
        for p in self.parameters: 
            #Check for missing required parameters: 
            if p.is_required and not(p.is_set): 
                raise ValueError("Parameter %s is not set." 
                                 % p.names[-1]) 
            #Also repeat the parameter validation here, just in case? 

    def _escape_filename(self, filename): 
        """Escape filenames with spaces by adding quotes (PRIVATE). 
        
        Note this will not add quotes if they are already included: 
        
        >>> print((_escape_filename('example with spaces'))) 
        "example with spaces" 
        >>> print((_escape_filename('"example with spaces"'))) 
        "example with spaces" 
        """ 
        # Is adding the following helpful 
        # if os.path.isfile(filename): 
        #    # On Windows, if the file exists, we can ask for 
        #    # its alternative short name (DOS style 8.3 format) 
        #    # which has no spaces in it.  Note that this name 
        #    # is not portable between machines, or even folder! 
        #    try: 
        #        import win32api 
        #        short = win32api.GetShortPathName(filename) 
        #        assert os.path.isfile(short) 
        #        return short 
        #    except ImportError: 
        #        pass 
        if " " not in filename: 
            return filename 
        # We'll just quote it - works on Windows, Mac OS X etc 
        if filename.startswith('"') and filename.endswith('"'): 
            # Its already quoted 
            return filename 
        else: 
            return '"%s"' % filename 

 
    def __str__(self): 
        """Make the commandline string with the currently set options. 
 
        e.g. 
        >>> from Bio.Emboss.Applications import WaterCommandline 
        >>> cline = WaterCommandline(gapopen=10, gapextend=0.5) 
        >>> cline.asequence = "asis:ACCCGGGCGCGGT" 
        >>> cline.bsequence = "asis:ACCCGAGCGCGGT" 
        >>> cline.outfile = "temp_water.txt" 
        >>> print(cline) 
        water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5 
        >>> str(cline) 
        'water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5' 
        """ 
        self._validate() 
        commandline = "%s " % self._escape_filename(self.program_name) 
        for parameter in self.parameters: 
           if parameter.is_set: 
                #This will include a trailing space: 
                commandline += str(parameter) 
        return commandline.strip() # remove trailing space 
 
    def __repr__(self): 
        """Return a representation of the command line object for debugging. 
 
        e.g. 
        >>> from Bio.Emboss.Applications import WaterCommandline 
        >>> cline = WaterCommandline(gapopen=10, gapextend=0.5) 
        >>> cline.asequence = "asis:ACCCGGGCGCGGT" 
        >>> cline.bsequence = "asis:ACCCGAGCGCGGT" 
        >>> cline.outfile = "temp_water.txt" 
        >>> print(cline) 
        water -outfile=temp_water.txt -asequence=asis:ACCCGGGCGCGGT -bsequence=asis:ACCCGAGCGCGGT -gapopen=10 -gapextend=0.5 
        >>> cline 
        WaterCommandline(cmd='water', outfile='temp_water.txt', asequence='asis:ACCCGGGCGCGGT', bsequence='asis:ACCCGAGCGCGGT', gapopen=10, gapextend=0.5) 
        """ 
        answer = "%s(cmd=%s" % (self.__class__.__name__, repr(self.program_name)) 
        for parameter in self.parameters: 
            if parameter.is_set: 
                if isinstance(parameter, _Switch): 
                    answer += ", %s=True" % parameter.names[-1] 
                else: 
                    answer += ", %s=%s" % (parameter.names[-1], repr(parameter.value)) 
        answer += ")" 
        return answer 
 
    def _get_parameter(self, name): 
        """Get a commandline option value.""" 
        for parameter in self.parameters: 
            if name in parameter.names: 
                if isinstance(parameter, _Switch): 
                    return parameter.is_set 
                else: 
                    return parameter.value 
        raise ValueError("Option name %s was not found." % name) 
 
    def _clear_parameter(self, name): 
        """Reset or clear a commandline option value.""" 
        cleared_option = False 
        for parameter in self.parameters: 
            if name in parameter.names: 
                parameter.value = None 
                parameter.is_set = False 
                cleared_option = True 
        if not cleared_option: 
            raise ValueError("Option name %s was not found." % name) 
 
    def set_parameter(self, name, value = None): 
        """Set a commandline option for a program (OBSOLETE). 
 
        Every parameter is available via a property and as a named 
        keyword when creating the instance. Using either of these is 
        preferred to this legacy set_parameter method which is now 
        OBSOLETE, and likely to be DEPRECATED and later REMOVED in 
        future releases. 
        """ 
        set_option = False 
        for parameter in self.parameters: 
            if name in parameter.names: 
                if isinstance(parameter, _Switch): 
                    if value is None: 
                        import warnings 
                        warnings.warn("For a switch type argument like %s, " 
                                      "we expect a boolean.  None is treated " 
                                      "as FALSE!" % parameter.names[-1]) 
                    parameter.is_set = bool(value) 
                    set_option = True 
                else: 
                    if value is not None: 
                        self._check_value(value, name, parameter.checker_function) 
                        parameter.value = value 
                    parameter.is_set = True 
                    set_option = True 
        if not set_option: 
            raise ValueError("Option name %s was not found." % name) 
 
    def _check_value(self, value, name, check_function): 
        """Check whether the given value is valid. 
 
        No return value - it either works or raises a ValueError. 
 
        This uses the passed function 'check_function', which can either 
        return a [0, 1] (bad, good) value or raise an error. Either way 
        this function will raise an error if the value is not valid, or 
        finish silently otherwise. 
        """ 
        if check_function is not None: 
            is_good = check_function(value)  # May raise an exception 
            assert is_good in [0, 1, True, False] 
            if not is_good: 
                raise ValueError("Invalid parameter value %r for parameter %s" 
                                 % (value, name)) 
 
    def __setattr__(self, name, value): 
        """Set attribute name to value (PRIVATE). 
 
        This code implements a workaround for a user interface issue. 
        Without this __setattr__ attribute-based assignment of parameters 
        will silently accept invalid parameters, leading to known instances 
        of the user assuming that parameters for the application are set, 
        when they are not. 
 
        >>> from Bio.Emboss.Applications import WaterCommandline 
        >>> cline = WaterCommandline(gapopen=10, gapextend=0.5, stdout=True) 
        >>> cline.asequence = "a.fasta" 
        >>> cline.bsequence = "b.fasta" 
        >>> cline.csequence = "c.fasta" 
        Traceback (most recent call last): 
        ... 
        ValueError: Option name csequence was not found. 
        >>> print(cline) 
        water -stdout -asequence=a.fasta -bsequence=b.fasta -gapopen=10 -gapextend=0.5 
 
        This workaround uses a whitelist of object attributes, and sets the 
        object attribute list as normal, for these.  Other attributes are 
        assumed to be parameters, and passed to the self.set_parameter method 
        for validation and assignment. 
        """ 
        if name in ['parameters', 'program_name']:  # Allowed attributes 
            self.__dict__[name] = value 
        else: 
            self.set_parameter(name, value)  # treat as a parameter 


 
    def __call__(self, stdin=None, stdout=True, stderr=True, 
                 cwd=None, env=None): 
        """Executes the command, waits for it to finish, and returns output. 
 
        Runs the command line tool and waits for it to finish. If it returns 
        a non-zero error level, an exception is raised. Otherwise two strings 
        are returned containing stdout and stderr. 
 
        The optional stdin argument should be a string of data which will be 
        passed to the tool as standard input. 
 
        The optional stdout and stderr argument may be filenames (string), 
        but otherwise are treated as a booleans, and control if the output 
        should be captured as strings (True, default), or ignored by sending 
        it to /dev/null to avoid wasting memory (False). If sent to a file 
        or ignored, then empty string(s) are returned. 
 
        The optional cwd argument is a string giving the working directory 
        to run the command from. See Python's subprocess module documentation 
        for more details. 
 
        The optional env argument is a dictionary setting the environment 
        variables to be used in the new process. By default the current 
        process' environment variables are used. See Python's subprocess 
        module documentation for more details. 
 
        Default example usage: 
 
        from Bio.Emboss.Applications import WaterCommandline 
        water_cmd = WaterCommandline(gapopen=10, gapextend=0.5, 
                                     stdout=True, auto=True, 
                                     asequence="a.fasta", bsequence="b.fasta") 
        print "About to run:\n%s" % water_cmd 
        std_output, err_output = water_cmd() 
 
        This functionality is similar to subprocess.check_output() added in 
        Python 2.7. In general if you require more control over running the 
        command, use subprocess directly. 
 
        As of Biopython 1.56, when the program called returns a non-zero error 
        level, a custom ApplicationError exception is raised. This includes 
        any stdout and stderr strings captured as attributes of the exception 
        object, since they may be useful for diagnosing what went wrong. 
        """ 
        if not stdout: 
            stdout_arg = open(os.devnull, "w") 
        elif isinstance(stdout, str): 
            stdout_arg = open(stdout, "w") 
        else: 
            stdout_arg = subprocess.PIPE 
 
        if not stderr: 
            stderr_arg = open(os.devnull, "w") 
        elif isinstance(stderr, str): 
            if stdout == stderr: 
                stderr_arg = stdout_arg #Write both to the same file 
            else: 
                stderr_arg = open(stderr, "w") 
        else: 
            stderr_arg = subprocess.PIPE 
 
        #We may not need to supply any piped input, but we setup the 
        #standard input pipe anyway as a work around for a python 
        #bug if this is called from a Windows GUI program.  For 
        #details, see http://bugs.python.org/issue1124861 
        # 
        #Using universal newlines is important on Python 3, this 
        #gives unicode handles rather than bytes handles. 
 
        #Windows 7 and 8 want shell = True 
        #platform is easier to understand that sys to determine 
        #windows version 
        if sys.platform != "win32": 
            use_shell = True 
        else: 
            win_ver = platform.win32_ver()[0] 
            if win_ver in ["7", "8"]: 
                use_shell = True 
            else: 
                use_shell = False 
        child_process = subprocess.Popen(str(self), stdin=subprocess.PIPE, stdout=stdout_arg, stderr=stderr_arg, universal_newlines=True, cwd=cwd, env=env, shell=use_shell) 

        #Use .communicate as can get deadlocks with .wait(), see Bug 2804 
        stdout_str, stderr_str = child_process.communicate(stdin) 
        if not stdout: 
            assert not stdout_str, stdout_str 
        if not stderr: 
            assert not stderr_str, stderr_str 
        return_code = child_process

        # Particularly important to close handles on Jython and PyPy 
        # (where garbage collection is less predictable) and on Windows 
        # (where cannot delete files with an open handle): 
        if not stdout or isinstance(stdout, str): 
            # We opened /dev/null or a file 
            stdout_arg.close() 
        if not stderr or (isinstance(stderr, str) and stdout != stderr): 
            # We opened /dev/null or a file 
            stderr_arg.close() 
        

        #I commented this out but I'm sure there should should be a reason for it. Fix it!
 #       if return_code: 
 #           raise ApplicationError(return_code, str(self), stdout_str, stderr_str)
        return stdout_str, stderr_str 



class _AbstractParameter: 
      """A class to hold information about a parameter for a commandline. 
      
      Do not use this directly, instead use one of the subclasses. 
      """ 
      def __init__(self): 
          raise NotImplementedError 
      
      def __str__(self): 
          raise NotImplementedError 
      
      
class _Option(_AbstractParameter): 
      """Represent an option that can be set for a program. 
      
      This holds UNIXish options like --append=yes and -a yes, 
      where a value (here "yes") is generally expected. 
      
      For UNIXish options like -kimura in clustalw which don't 
      take a value, use the _Switch object instead. 
      
      Attributes: 
      
      o names -- a list of string names (typically two entries) by which 
      the parameter can be set via the legacy set_parameter method 
      (eg ["-a", "--append", "append"]). The first name in list is used 
      when building the command line. The last name in the list is a 
      "human readable" name describing the option in one word. This 
      must be a valid Python identifier as it is used as the property 
      name and as a keyword argument, and should therefore follow PEP8 
      
      o description -- a description of the option. This is used as 
      the property docstring. 
      
      o filename -- True if this argument is a filename and should be 
      automatically quoted if it contains spaces. 
      
      o checker_function -- a reference to a function that will determine 
      if a given value is valid for this parameter. This function can either 
      raise an error when given a bad value, or return a [0, 1] decision on 
      whether the value is correct. 
      
      o equate -- should an equals sign be inserted if a value is used? 
      
      o is_required -- a flag to indicate if the parameter must be set for 
      the program to be run. 
      
      o is_set -- if the parameter has been set 
      
      o value -- the value of a parameter 
      """ 
      def __init__(self, names, description, filename=False, checker_function=None, 
                   is_required=False, equate=True): 
          self.names = names 
          assert isinstance(description, str), "%r for %s" % (description, names[-1]) 
          self.is_filename = filename 
          self.checker_function = checker_function 
          self.description = description 
          self.equate = equate 
          self.is_required = is_required 
      
          self.is_set = False 
          self.value = None 
      
      def __str__(self): 
          """Return the value of this option for the commandline. 
      
          Includes a trailing space. 
          """ 
          # Note: Before equate was handled explicitly, the old 
          # code would do either "--name " or "--name=value ", 
          # or " -name " or " -name value ".  This choice is now 
          # now made explicitly when setting up the option. 
          if self.value is None: 
              return "%s " % self.names[0] 
          if self.is_filename: 
              v = _escape_filename(self.value) 
          else: 
              v = str(self.value) 
          if self.equate: 
              return "%s=%s " % (self.names[0], v) 
          else: 
              return "%s %s " % (self.names[0], v) 
      
      
class _Switch(_AbstractParameter): 
      """Represent an optional argument switch for a program. 
      
      This holds UNIXish options like -kimura in clustalw which don't 
      take a value, they are either included in the command string 
      or omitted. 
      
      o names -- a list of string names (typically two entries) by which 
      the parameter can be set via the legacy set_parameter method 
      (eg ["-a", "--append", "append"]). The first name in list is used 
      when building the command line. The last name in the list is a 
      "human readable" name describing the option in one word. This 
      must be a valid Python identifer as it is used as the property 
      name and as a keyword argument, and should therefore follow PEP8 
      naming. 
      
      o description -- a description of the option. This is used as 
      the property docstring. 
      
      o is_set -- if the parameter has been set 
      
      NOTE - There is no value attribute, see is_set instead, 
      """ 
      def __init__(self, names, description): 
          self.names = names 
          self.description = description 
          self.is_set = False 
          self.is_required = False 
      
      def __str__(self): 
          """Return the value of this option for the commandline. 
      
          Includes a trailing space. 
          """ 
          assert not hasattr(self, "value") 
          if self.is_set: 
              return "%s " % self.names[0] 
          else: 
              return "" 
      
      
class _Argument(_AbstractParameter): 
      """Represent an argument on a commandline. 
      
      The names argument should be a list containing one string. 
      This must be a valid Python identifer as it is used as the 
      property name and as a keyword argument, and should therefore 
      follow PEP8 naming. 
      """ 
      def __init__(self, names, description, filename=False, 
                   checker_function=None, is_required=False): 
          # if len(names) != 1: 
          #    raise ValueError("The names argument to _Argument should be a " 
          #                     "single entry list with a PEP8 property name.") 
          self.names = names 
          assert isinstance(description, str), "%r for %s" % (description, names[-1]) 
          self.is_filename = filename 
          self.checker_function = checker_function 
          self.description = description 
          self.is_required = is_required 
          self.is_set = False 
          self.value = None 
         
      def __str__(self): 
          if self.value is None: 
              return " " 
          elif self.is_filename: 
              return "%s " % _escape_filename(self.value) 
          else: 
              return "%s " % self.value 
      
      
class _ArgumentList(_Argument): 
      """Represent a variable list of arguments on a command line, e.g. multiple filenames.""" 
      # TODO - Option to require at least one value? e.g. min/max count? 
      
      def __str__(self): 
          assert isinstance(self.value, list), "Arguments should be a list" 
          assert self.value, "Requires at least one filename" 
          # A trailing space is required so that parameters following the last filename 
          # do not appear merged. 
          # e.g.:  samtools cat in1.bam in2.bam-o out.sam  [without trailing space][Incorrect] 
          #        samtools cat in1.bam in2.bam -o out.sam  [with trailing space][Correct] 
          if self.is_filename: 
              return " ".join(_escape_filename(v) for v in self.value) + " " 
          else: 
              return " ".join(self.value) + " " 
      
      
class _StaticArgument(_AbstractParameter): 
      """Represent a static (read only) argument on a commandline. 
      
      This is not intended to be exposed as a named argument or 
      property of a command line wrapper object. 
      """ 
      def __init__(self, value): 
          self.names = [] 
          self.is_required = False 
          self.is_set = True 
          self.value = value 
      
      def __str__(self): 
          return "%s " % self.value 
      
      
def _escape_filename(filename): 
      """Escape filenames with spaces by adding quotes (PRIVATE). 
      
      Note this will not add quotes if they are already included: 
      
      >>> print((_escape_filename('example with spaces'))) 
      "example with spaces" 
      >>> print((_escape_filename('"example with spaces"'))) 
      "example with spaces" 
      """ 
      # Is adding the following helpful 
      # if os.path.isfile(filename): 
      #    # On Windows, if the file exists, we can ask for 
      #    # its alternative short name (DOS style 8.3 format) 
      #    # which has no spaces in it.  Note that this name 
      #    # is not portable between machines, or even folder! 
      #    try: 
      #        import win32api 
      #        short = win32api.GetShortPathName(filename) 
      #        assert os.path.isfile(short) 
      #        return short 
      #    except ImportError: 
      #        pass 
      if " " not in filename: 
          return filename 
      # We'll just quote it - works on Windows, Mac OS X etc 
      if filename.startswith('"') and filename.endswith('"'): 
          # Its already quoted 
          return filename 
      else: 
          return '"%s"' % filename 











class MuscleCommandline(AbstractCommandline): 
    r"""Command line wrapper for the multiple alignment program MUSCLE. 
 
    http://www.drive5.com/muscle/ 
 
    Example: 
    -------- 
 
    >>> from Bio.Align.Applications import MuscleCommandline 
    >>> muscle_exe = r"C:\Program Files\Aligments\muscle3.8.31_i86win32.exe" 
    >>> in_file = r"C:\My Documents\unaligned.fasta" 
    >>> out_file = r"C:\My Documents\aligned.fasta" 
    >>> muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file) 
    >>> print(muscle_cline) 
    "C:\Program Files\Aligments\muscle3.8.31_i86win32.exe" -in "C:\My Documents\unaligned.fasta" -out "C:\My Documents\aligned.fasta" 
 
    You would typically run the command line with muscle_cline() or via 
    the Python subprocess module, as described in the Biopython tutorial. 
 
    Citations: 
    ---------- 
 
    Edgar, Robert C. (2004), MUSCLE: multiple sequence alignment with high 
    accuracy and high throughput, Nucleic Acids Research 32(5), 1792-97. 
 
    Edgar, R.C. (2004) MUSCLE: a multiple sequence alignment method with 
    reduced time and space complexity. BMC Bioinformatics 5(1): 113. 
 
    Last checked against version: 3.7, briefly against 3.8 
    """ 
    def __init__(self, cmd="muscle", **kwargs): 
        CLUSTERING_ALGORITHMS = ["upgma", "upgmb", "neighborjoining"] 
        DISTANCE_MEASURES_ITER1 = ["kmer6_6", "kmer20_3", "kmer20_4", "kbit20_3", 
                                   "kmer4_6"] 
        DISTANCE_MEASURES_ITER2 = DISTANCE_MEASURES_ITER1 + ["pctid_kimura", "pctid_log"] 
        OBJECTIVE_SCORES = ["sp", "ps", "dp", "xp", "spf", "spm"] 
        TREE_ROOT_METHODS = ["pseudo", "midlongestspan", "minavgleafdist"] 
        SEQUENCE_TYPES = ["protein", "nucleo", "auto"] 
        WEIGHTING_SCHEMES = ["none", "clustalw", "henikoff", "henikoffpb", "gsc", "threeway"] 
        self.parameters = [ 
            # Can't use "in" as the final alias as this is a reserved word in python: 
            _Option(["-in", "in", "input"], 
                    "Input filename", 
                    filename=True, 
                    equate=False), 
            _Option(["-out", "out"], 
                    "Output filename", 
                    filename=True, 
                    equate=False), 
            _Switch(["-diags", "diags"], 
                    "Find diagonals (faster for similar sequences)"), 
            _Switch(["-profile", "profile"], 
                    "Perform a profile alignment"), 
            _Option(["-in1", "in1"], 
                    "First input filename for profile alignment", 
                    filename=True, 
                    equate=False), 
            _Option(["-in2", "in2"], 
                    "Second input filename for a profile alignment", 
                    filename=True, 
                    equate=False), 
            # anchorspacing   Integer              32                 Minimum spacing between 
            _Option(["-anchorspacing", "anchorspacing"], 
                    "Minimum spacing between anchor columns", 
                    checker_function=lambda x: isinstance(x, int), 
                    equate=False), 
            # center          Floating point       [1]                Center parameter. 
            #                                                        Should be negative. 
            _Option(["-center", "center"], 
                    "Center parameter - should be negative", 
                    checker_function=lambda x: isinstance(x, float), 
                    equate=False), 
            # cluster1        upgma                upgmb              Clustering method. 
            _Option(["-cluster1", "cluster1"], 
                    "Clustering method used in iteration 1", 
                    checker_function=lambda x: x in CLUSTERING_ALGORITHMS, 
                    equate=False), 
            # cluster2        upgmb                                   cluster1 is used in 
            #                neighborjoining                         iteration 1 and 2, 
            #                                                        cluster2 in later 
            #                                                        iterations. 
            _Option(["-cluster2", "cluster2"], 
                    "Clustering method used in iteration 2", 
                    checker_function=lambda x: x in CLUSTERING_ALGORITHMS, 
                    equate=False), 
            # diaglength      Integer              24                 Minimum length of 
            #                                                        diagonal. 
            _Option(["-diaglength", "diaglength"], 
                    "Minimum length of diagonal", 
                    checker_function=lambda x: isinstance(x, int), 
                    equate=True), 
            # diagmargin      Integer              5                  Discard this many 
            #                                                        positions at ends of 
            #                                                        diagonal. 
            _Option(["-diagmargin", "diagmargin"], 
                    "Discard this many positions at ends of diagonal", 
                    checker_function=lambda x: isinstance(x, int), 
                    equate=False), 
            # distance1       kmer6_6              Kmer6_6 (amino) or Distance measure for 
            #                kmer20_3             Kmer4_6 (nucleo)   iteration 1. 
            #                kmer20_4 
            #                kbit20_3 
            #                kmer4_6 
            _Option(["-distance1", "distance1"], 
                    "Distance measure for iteration 1", 
                    checker_function=lambda x: x in DISTANCE_MEASURES_ITER1, 
                    equate=False), 
            # distance2       kmer6_6              pctid_kimura       Distance measure for 
            #                kmer20_3                                iterations 2, 3 ... 
            #                kmer20_4 
            #                kbit20_3 
            #                pctid_kimura 
            #                pctid_log 
            _Option(["-distance2", "distance2"], 
                    "Distance measure for iteration 2", 
                    checker_function=lambda x: x in DISTANCE_MEASURES_ITER2, 
                    equate=False), 
            # gapopen         Floating point       [1]                The gap open score. 
            #                                                        Must be negative. 
            _Option(["-gapopen", "gapopen"], 
                    "Gap open score - negative number", 
                    checker_function=lambda x: isinstance(x, float), 
                    equate=False), 
            # hydro           Integer              5                  Window size for 
            #                                                        determining whether a 
            #                                                        region is hydrophobic. 
            _Option(["-hydro", "hydro"], 
                    "Window size for hydrophobic region", 
                    checker_function=lambda x: isinstance(x, int), 
                    equate=False), 
            # hydrofactor     Floating point       1.2                Multiplier for gap 
            #                                                        open/close penalties in 
            #                                                        hydrophobic regions. 
            _Option(["-hydrofactor", "hydrofactor"], 
                    "Multiplier for gap penalties in hydrophobic regions", 
                    checker_function=lambda x: isinstance(x, float), 
                    equate=False), 
            # log             File name            None.              Log file name (delete 
            #                                                        existing file). 
            _Option(["-log", "log"], 
                    "Log file name", 
                    filename=True, 
                    equate=False), 
            # loga            File name            None.              Log file name (append 
            #                                                        to existing file). 
            _Option(["-loga", "loga"], 
                    "Log file name (append to existing file)", 
                    filename=True, 
                    equate=False), 
            # maxdiagbreak    Integer              1                  Maximum distance 
            #                                                        between two diagonals 
            #                                                        that allows them to 
            #                                                        merge into one 
            #                                                        diagonal. 
            _Option(["-maxdiagbreak", "maxdiagbreak"], 
                    "Maximum distance between two diagonals that allows " 
                    "them to merge into one diagonal", 
                    checker_function=lambda x: isinstance(x, int), 
                    equate=False), 
            # maxhours        Floating point       None.              Maximum time to run in 
            #                                                        hours. The actual time 
            #                                                        may exceed the 
            #                                                        requested limit by a 
            #                                                        few minutes. Decimals 
            #                                                        are allowed, so 1.5 
            #                                                        means one hour and 30 
            #                                                        minutes. 
            _Option(["-maxhours", "maxhours"], 
                    "Maximum time to run in hours", 
                    checker_function=lambda x: isinstance(x, float), 
                    equate=False), 
            # maxiters        Integer 1, 2 ...     16                 Maximum number of 
            #                                                        iterations. 
            _Option(["-maxiters", "maxiters"], 
                    "Maximum number of iterations", 
                    checker_function=lambda x: isinstance(x, int), 
                    equate=False), 
            # maxtrees        Integer              1                  Maximum number of new 
            #                                                        trees to build in 
            #                                                        iteration 2. 
            _Option(["-maxtrees", "maxtrees"], 
                    "Maximum number of trees to build in iteration 2", 
                    checker_function=lambda x: isinstance(x, int), 
                    equate=False), 
            # minbestcolscore Floating point       [1]                Minimum score a column 
            #                                                        must have to be an 
            #                                                        anchor. 
            _Option(["-minbestcolscore", "minbestcolscore"], 
                    "Minimum score a column must have to be an anchor", 
                    checker_function=lambda x: isinstance(x, float), 
                    equate=False), 
            # minsmoothscore  Floating point       [1]                Minimum smoothed score 
            #                                                        a column must have to 
            #                                                        be an anchor. 
            _Option(["-minsmoothscore", "minsmoothscore"], 
                    "Minimum smoothed score a column must have to " 
                    "be an anchor", 
                    checker_function=lambda x: isinstance(x, float), 
                    equate=False), 
            # objscore        sp                   spm                Objective score used by 
            #                ps                                      tree dependent 
            #                dp                                      refinement. 
            #                xp                                      sp=sum-of-pairs score. 
            #                spf                                     spf=sum-of-pairs score 
            #                spm                                     (dimer approximation) 
            #                                                        spm=sp for < 100 seqs, 
            #                                                        otherwise spf 
            #                                                        dp=dynamic programming 
            #                                                        score. 
            #                                                        ps=average profile- 
            #                                                        sequence score. 
            #                                                        xp=cross profile score. 
            _Option(["-objscore", "objscore"], 
                    "Objective score used by tree dependent refinement", 
                    checker_function=lambda x: x in OBJECTIVE_SCORES, 
                    equate=False), 
            # root1           pseudo               pseudo             Method used to root 
            _Option(["-root1", "root1"], 
                    "Method used to root tree in iteration 1", 
                    checker_function=lambda x: x in TREE_ROOT_METHODS, 
                    equate=False), 
            # root2           midlongestspan                          tree; root1 is used in 
            #                minavgleafdist                          iteration 1 and 2, 
            #                                                        root2 in later 
            #                                                        iterations. 
            _Option(["-root2", "root2"], 
                    "Method used to root tree in iteration 2", 
                    checker_function=lambda x: x in TREE_ROOT_METHODS, 
                    equate=False), 
            # seqtype         protein              auto               Sequence type. 
            #                nucleo 
            #                auto 
            _Option(["-seqtype", "seqtype"], 
                    "Sequence type", 
                    checker_function=lambda x: x in SEQUENCE_TYPES, 
                    equate=False), 
            # smoothscoreceil Floating point       [1]                Maximum value of column 
            #                                                        score for smoothing 
            #                                                        purposes. 
            _Option(["-smoothscoreceil", "smoothscoreceil"], 
                    "Maximum value of column score for smoothing", 
                    checker_function=lambda x: isinstance(x, float), 
                    equate=False), 
            # smoothwindow    Integer              7                  Window used for anchor 
            #                                                        column smoothing. 
            _Option(["-smoothwindow", "smoothwindow"], 
                    "Window used for anchor column smoothing", 
                    checker_function=lambda x: isinstance(x, int), 
                    equate=False), 
            # SUEFF           Floating point value 0.1                Constant used in UPGMB 
            #                between 0 and 1.                        clustering. Determines 
            #                                                        the relative fraction 
            #                                                        of average linkage 
            #                                                        (SUEFF) vs. nearest- 
            #                                                        neighbor linkage (1 
            #                                                        SUEFF). 
            _Option(["-sueff", "sueff"], 
                    "Constant used in UPGMB clustering", 
                    checker_function=lambda x: isinstance(x, float), 
                    equate=False), 
            # tree1           File name            None               Save tree produced in 
            _Option(["-tree1", "tree1"], 
                    "Save Newick tree from iteration 1", 
                    equate=False), 
            # tree2                                                   first or second 
            #                                                        iteration to given file 
            #                                                        in Newick (Phylip- 
            #                                                        compatible) format. 
            _Option(["-tree2", "tree2"], 
                    "Save Newick tree from iteration 2", 
                    equate=False), 
            # weight1         none                 clustalw           Sequence weighting 
            _Option(["-weight1", "weight1"], 
                    "Weighting scheme used in iteration 1", 
                    checker_function=lambda x: x in WEIGHTING_SCHEMES, 
                    equate=False), 
            # weight2         henikoff                                scheme. 
            #                henikoffpb                              weight1 is used in 
            #                gsc                                     iterations 1 and 2. 
            #                clustalw                                weight2 is used for 
            #                threeway                                tree-dependent 
            #                                                        refinement. 
            #                                                        none=all sequences have 
            #                                                        equal weight. 
            #                                                        henikoff=Henikoff & 
            #                                                        Henikoff weighting 
            #                                                        scheme. 
            #                                                        henikoffpb=Modified 
            #                                                        Henikoff scheme as used 
            #                                                        in PSI-BLAST. 
            #                                                        clustalw=CLUSTALW 
            #                                                        method. 
            #                                                        threeway=Gotoh three- 
            #                                                        way method. 
            _Option(["-weight2", "weight2"], 
                    "Weighting scheme used in iteration 2", 
                    checker_function=lambda x: x in WEIGHTING_SCHEMES, 
                    equate=False), 
            # ################### FORMATS ####################################### 
            # Multiple formats can be specified on the command line 
            # If -msf appears it will be used regardless of other formats 
            # specified. If -clw appears (and not -msf), clustalw format will be 
            # used regardless of other formats specified. If both -clw and 
            # -clwstrict are specified -clwstrict will be used regardless of 
            # other formats specified. If -fasta is specified and not -msf, 
            # -clw, or clwstrict, fasta will be used. If -fasta and -html are 
            # specified -fasta will be used. Only if -html is specified alone 
            # will html be used. I kid ye not. 
            # clw                no              Write output in CLUSTALW format (default is 
            #                                   FASTA). 
            _Switch(["-clw", "clw"], 
                    "Write output in CLUSTALW format (with a MUSCLE header)"), 
            # clwstrict          no              Write output in CLUSTALW format with the 
            #                                   "CLUSTAL W (1.81)" header rather than the 
            #                                   MUSCLE version. This is useful when a post- 
            #                                   processing step is picky about the file 
            #                                   header. 
            _Switch(["-clwstrict", "clwstrict"], 
                    "Write output in CLUSTALW format with version 1.81 header"), 
            # fasta              yes             Write output in FASTA format. Alternatives 
            #                                   include clw, 
            #                                   clwstrict, msf and html. 
            _Switch(["-fasta", "fasta"], 
                    "Write output in FASTA format"), 
            # html               no              Write output in HTML format (default is 
            #                                   FASTA). 
            _Switch(["-html", "html"], 
                    "Write output in HTML format"), 
            # msf                no              Write output in MSF format (default is 
            #                                   FASTA). 
            _Switch(["-msf", "msf"], 
                    "Write output in MSF format"), 
            # Phylip interleaved - undocumented as of 3.7 
            _Switch(["-phyi", "phyi"], 
                    "Write output in PHYLIP interleaved format"), 
            # Phylip sequential - undocumented as of 3.7 
            _Switch(["-phys", "phys"], 
                    "Write output in PHYLIP sequential format"), 
            # ################# Additional specified output files ######### 
            _Option(["-phyiout", "phyiout"], 
                    "Write PHYLIP interleaved output to specified filename", 
                    filename=True, 
                    equate=False), 
            _Option(["-physout", "physout"], "Write PHYLIP sequential format to specified filename", 
                    filename=True, 
                    equate=False), 
            _Option(["-htmlout", "htmlout"], "Write HTML output to specified filename", 
                    filename=True, 
                    equate=False), 
            _Option(["-clwout", "clwout"], 
                    "Write CLUSTALW output (with MUSCLE header) to specified " 
                    "filename", 
                    filename=True, 
                    equate=False), 
            _Option(["-clwstrictout", "clwstrictout"], 
                    "Write CLUSTALW output (with version 1.81 header) to " 
                    "specified filename", 
                    filename=True, 
                    equate=False), 
            _Option(["-msfout", "msfout"], 
                    "Write MSF format output to specified filename", 
                    filename=True, 
                    equate=False), 
            _Option(["-fastaout", "fastaout"], 
                    "Write FASTA format output to specified filename", 
                    filename=True, 
                    equate=False), 
            # ############# END FORMATS ################################### 
            # anchors            yes             Use anchor optimization in tree dependent 
            #                                   refinement iterations. 
            _Switch(["-anchors", "anchors"], 
                    "Use anchor optimisation in tree dependent " 
                    "refinement iterations"), 
            # noanchors          no              Disable anchor optimization. Default is 
            #                                   anchors. 
            _Switch(["-noanchors", "noanchors"], 
                    "Do not use anchor optimisation in tree dependent " 
                    "refinement iterations"), 
            # group              yes             Group similar sequences together in the 
            #                                   output. This is the default. See also 
            #                                   stable. 
            _Switch(["-group", "group"], 
                    "Group similar sequences in output"), 
            # stable             no              Preserve input order of sequences in output 
            #                                   file. Default is to group sequences by 
            #                                   similarity (group). 
            _Switch(["-stable", "stable"], 
                    "Do not group similar sequences in output (not supported in v3.8)"), 
            # ############# log-expectation profile score ###################### 
            # One of either -le, -sp, or -sv 
            # 
            # According to the doc, spn is default and the only option for 
            # nucleotides: this doesnt appear to be true. -le, -sp, and -sv can 
            # be used and produce numerically different logs (what is going on?) 
            # 
            # spn fails on proteins 
            # le                 maybe           Use log-expectation profile score (VTML240). 
            #                                    Alternatives are to use sp or sv. This is 
            #                                    the default for amino acid sequences. 
            _Switch(["-le", "le"], 
                    "Use log-expectation profile score (VTML240)"), 
            # sv                 no              Use sum-of-pairs profile score (VTML240). 
            #                                   Default is le. 
            _Switch(["-sv", "sv"], 
                    "Use sum-of-pairs profile score (VTML240)"), 
            # sp                 no              Use sum-of-pairs protein profile score 
            #                                   (PAM200). Default is le. 
            _Switch(["-sp", "sp"], 
                    "Use sum-of-pairs protein profile score (PAM200)"), 
            # spn                maybe           Use sum-of-pairs nucleotide profile score 
            #                                   (BLASTZ parameters). This is the only option 
            #                                   for nucleotides, and is therefore the 
            #                                   default. 
            _Switch(["-spn", "spn"], 
                    "Use sum-of-pairs protein nucleotide profile score"), 
            # ############# END log-expectation profile score ###################### 
            # quiet              no              Do not display progress messages. 
            _Switch(["-quiet", "quiet"], 
                    "Use sum-of-pairs protein nucleotide profile score"), 
            # refine             no              Input file is already aligned, skip first 
            #                                   two iterations and begin tree dependent 
            #                                   refinement. 
            _Switch(["-refine", "refine"], 
                    "Only do tree dependent refinement"), 
            # core               yes in muscle,  Do not catch exceptions. 
            #                   no in muscled. 
            _Switch(["-core", "core"], 
                    "Catch exceptions"), 
            # nocore             no in muscle,   Catch exceptions and give an error message 
            #                   yes in muscled. if possible. 
            _Switch(["-nocore", "nocore"], 
                    "Do not catch exceptions"), 
            # termgapsfull       no              Terminal gaps penalized with full penalty. 
            #                                   [1] Not fully supported in this version. 
            # 
            # termgapshalf       yes             Terminal gaps penalized with half penalty. 
            #                                   [1] Not fully supported in this version. 
            # 
            # termgapshalflonger no              Terminal gaps penalized with half penalty if 
            #                                   gap relative to 
            #                                   longer sequence, otherwise with full 
            #                                   penalty. 
            #                                   [1] Not fully supported in this version. 
            # verbose            no              Write parameter settings and progress 
            #                                   messages to log file. 
            _Switch(["-verbose", "verbose"], 
                    "Write parameter settings and progress"), 
            # version            no              Write version string to stdout and exit. 
            _Switch(["-version", "version"], 
                    "Write version string to stdout and exit"), 
           ] 
        AbstractCommandline.__init__(self, cmd, **kwargs) 
 
 

