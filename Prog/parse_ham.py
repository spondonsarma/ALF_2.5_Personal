#!/usr/bin/env python3
"""Script for automatically parsing parammeters of hamiltonian"""
# pylint: disable=invalid-name

# import os
import sys
import pprint


def parse(filename):
    """Parse Fortran file for parameter lists.

    Each list has to start with a line containing '#PARAMETERS START#' and end
    with a line containing '#PARAMETERS END#'. The name written in the same
    line as '#PARAMETERS START#' is the namelinst name."""
    parameters = {}
    with open(filename, 'r') as f:
        lines = f.readlines()

    do_parse = False
    for line in lines:
        if '#PARAMETERS START#' in line:
            do_parse = True
            name_key = line.split('#PARAMETERS START#')[1].strip()
            namelist = {}
            parameters[name_key] = namelist
            continue
        if '#PARAMETERS END#' in line:
            do_parse = False
            continue
        if do_parse and '::' in line:
            par = parse_line(line)
            par_name = par.pop('name')
            namelist[par_name] = par
    return parameters


def parse_line(line):
    """Parse single line in Fortran file for parameter."""
    parameter = {}
    dtype, rest = line.split('::', maxsplit=1)
    dtype = dtype.strip().lower()
    if '!' in rest:
        assignment, comment = rest.strip().split('!', maxsplit=1)
        assignment = assignment.strip()
        parameter['comment'] = comment.strip(' !')
    else:
        assignment = rest.strip()
        parameter['comment'] = ''
    name, value = assignment.split('=')
    # parameter['name'] = name.split('(')[0].strip()
    parameter['name'] = name.strip()
    value = value.strip().replace('d0', '0')
    parameter['defined_in_base'] = dtype.startswith('!')
    if '[' in value:
        parameter['value'] = []
        for val in value.strip('[]').split(','):
            parameter['value'].append(_to_value(dtype, val))
    else:
        parameter['value'] = _to_value(dtype, value)

    return parameter


def _max_len(dictionary):
    return max([len(i) for i in dictionary])


def _dtype_name(parameter):
    if isinstance(parameter, list):
        return _dtype_name(parameter[0])
    if isinstance(parameter, bool):
        return 'logical'
    if isinstance(parameter, float):
        return 'real(dp)'
    if isinstance(parameter, int):
        return 'integer'
    if isinstance(parameter, str):
        return 'Character(len=64)'
    raise Exception('Error in "_dtype_name": unrecognized type')


def _convert_par_to_str(parameter):
    """Converts a given parameter value to a string that can be
    written into a parameter file.
    """
    if isinstance(parameter, list):
        pars = [_convert_par_to_str(p) for p in parameter]
        return '[{}]'.format(', '.join(pars))
    if isinstance(parameter, bool):
        if parameter:
            return '.true.'
        return '.false.'
    if isinstance(parameter, float):
        if 'e' in '{}'.format(parameter):
            return '{}'.format(parameter).replace('e', 'd')
        return '{}d0'.format(parameter)
    if isinstance(parameter, int):
        return '{}'.format(parameter)
    if isinstance(parameter, str):
        return '"{}"'.format(parameter)

    raise Exception('Error in "_convert_par_to_str": unrecognized type')


def create_read_par(filename, parameters):
    INDENT = 9     # Number of indentation spaces
    LINE_MAX = 70  # Maximal line length
    with open('read_parameters_template', 'r') as f:
        lines = f.readlines()

    f = open(filename, 'w')
    for line in lines:
        if '##NAMELIST##' in line:
            for nlist_name, nlist in parameters.items():
                s = '{}NAMELIST /{}/  '.format(INDENT*' ', nlist_name)
                for par_name in nlist:
                    par_name = par_name.split('(')[0]
                    if len(s) < LINE_MAX:
                        s = '{}{}, '.format(s, par_name)
                    else:
                        f.write('{}&\n'.format(s))
                        s = '{}     &     {}, '.format(INDENT*' ', par_name)
                f.write(s[:-2]+'\n')
        elif '##PARAMETER_DEF##' in line:
            for nlist_name, nlist in parameters.items():
                f.write('{}!Parameters {}\n'.format(INDENT*' ', nlist_name))

                names_len = _max_len(nlist)
                # dtypes_str = [_dtype_name(par['value']) for
                #               par_name, par in nlist.items()]
                # dtypes_len = _max_len(dtypes_str)
                pars_str = [_convert_par_to_str(par['value']) for
                            par_name, par in nlist.items()]
                pars_len = _max_len(pars_str)
                # comments = [par[1] for par_name, par in nlist.items()]
                # for i in range(len(nlist)):

                for par_name, par in nlist.items():
                    # s = '{}{} :: {} = {} !'.format(
                    #     INDENT*' ',
                    #     _dtype_name(par['value']).ljust(dtypes_len),
                    #     par_name.ljust(names_len),
                    #     _convert_par_to_str(par['value']).ljust(pars_len)
                    #     )
                    s = '{}{} = {} !'.format(
                        INDENT*' ',
                        par_name.ljust(names_len),
                        _convert_par_to_str(par['value']).ljust(pars_len)
                        )
                    comment = par['comment'].split(' ')
                    comment_indent = min(40, len(s)-2)
                    for word in comment:
                        if len(s) < LINE_MAX:
                            s = '{} {}'.format(s, word)
                        else:
                            f.write('{} \n'.format(s))
                            s = '{} ! {}'.format(comment_indent*' ', word)
                    f.write(s)
                    f.write('\n')
        elif '##READ_VAR##' in line:
            for nlist_name, nlist in parameters.items():
                s = '{}   READ(unit_para, NML={})\n'.format(INDENT*' ', nlist_name)
                f.write(s)

        # Relevant if parameters in hamiltonian type
        # elif '##COPY_IN_HAM##' in line:
        #     for nlist_name, nlist in parameters.items():
        #         for par_name in nlist:
        #             s = '{0:}   ham%{1:} = {1:}\n'.format(
        #                 INDENT*' ',
        #                 par_name.split('(')[0].ljust(names_len))
        #             f.write(s)

        elif '##MPI_BCAST##' in line:
            fstring = '{}CALL MPI_BCAST({},{:>3},{},0,Group_Comm,ierr)\n'
            for nlist_name, nlist in parameters.items():
                names_len = _max_len(nlist)
                for par_name, par in nlist.items():
                    s = fstring.format(
                        INDENT*' ',
                        par_name.split('(')[0].ljust(names_len),
                        _get_mpi_len(par['value']),
                        _get_mpi_dtype(par['value']),
                        )
                    f.write(s)
        else:
            f.write(line)
    f.close()


def _to_value(dtype, value):
    if 'real' in dtype:
        return float(value)
    if 'integer' in dtype:
        return int(value)
    if 'character' in dtype:
        return value.strip('"\'')
    if 'logical' in dtype:
        if 't' in value or 'T' in value:
            return True
        if 'f' in value or 'F' in value:
            return False

        raise Exception(
            '"{}" can not be mapped to bool.'.format(value))

    raise Exception('"{}" can not be mapped to a type'.format(dtype))


def _get_mpi_dtype(parameter):
    if isinstance(parameter, list):
        return _get_mpi_dtype(parameter[0])
    if isinstance(parameter, bool):
        return 'MPI_LOGICAL  '
    if isinstance(parameter, float):
        return 'MPI_REAL8    '
    if isinstance(parameter, int):
        return 'MPI_INTEGER  '
    if isinstance(parameter, str):
        return 'MPI_CHARACTER'
    raise Exception('Error in "_get_mpi_dtype": unrecognized type')


def _get_mpi_len(parameter):
    if isinstance(parameter, list):
        return _get_mpi_len(parameter[0]) * len(parameter)
    if isinstance(parameter, (bool, float, int)):
        return 1
    if isinstance(parameter, str):
        return 64
    raise Exception('Error in "_get_mpi_dtype": unrecognized type')


def parse_list(ham_names_file):
    """."""
    parameters = {}
    with open(ham_names_file, 'r') as f:
        ham_names = f.read().splitlines()

    for ham_name in ham_names:
        filename = 'Hamiltonians/Hamiltonian_{}_smod.F90'.format(ham_name)
        print('filename:', filename)

        parameters = parse(filename)
        pprint.pprint(parameters)
        filename = 'Hamiltonians/Hamiltonian_{}_read_parameters.F90'.format(ham_name)
        create_read_par(filename, parameters)


if __name__ == '__main__':
    parse_list(sys.argv[1])
