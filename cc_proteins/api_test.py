import requests
import json

API_UNITS = 'http://127.0.0.1:9000/api/units'
API_CONVERT = 'http://127.0.0.1:9000/api/convert'

# simple function to stringify a list
def iterable_to_str(xs, sep=' '):
    return sep.join([str(x) for x in xs])


def call_api(url, data=None):
    resp = requests.post(url,
                         data=json.dumps(data),
                         headers={'Content-type': 'application/json',
                                  'Accept': 'application/json'})

    if resp.status_code == 200:
        json_data = json.loads(resp.text)
        return json_data
    else:
        raise ValueError(resp.text)

# do a conversion using the JSON rest service and print out the results
def convert(req_dict):
    # prepare and POST the JSON request
    #print('converting from "{}" to "{}":'.format(req_dict['fromUnit'], req_dict['toUnit']))
    json_data = call_api(API_CONVERT, req_dict)
    return json_data['results']


def get_units():
    units = call_api(API_UNITS)
    return units['conversion']


def main():
    units = get_units()

    bp_coordinates = [
        ['1', 92887845],
        ['1', 172956646],
        ['X', 8574635]
    ]
    '''
    ['1', 17299956646],
    ['2', 13126731],
    ['3', 7821365],
    ['4', 41364422],
    ['5', 92572826],
    ['6', 47433029],
    ['Chr7', 61246731],
    ['8', 31246731],
    ['9', 61246731],
    ['10', 60246731],
    ['11', 91246731],
    ['12', 61236831],
    ['13', 60246291],
    ['14', 41246181],
    ['15', 20446772],
    ['16', 61245863],
    ['17', 82244954],
    ['18', 19243145],
    ['19', 70242036],
    '''

    cm_coordinates = [
        ['1', 34.2],
        ['19', 20.24],
        ['X', 14.2]
    ]
    '''
        ['1', 0],
        ['1', 87.3],
        ['2', 13126731.1],
        ['2', 45.2],
        ['3', 11.6],
        ['4', 12.982],
        ['5', 34.7],
        ['6', 44.89],
        ['Chr7', 62.3],
        ['8', 48.0123],
        ['9', 12.345],
        ['10', 67.890],
        ['11', 11.72],
        ['12', 68.3423],
        ['13', 27.23],
        ['14', 11.11],
        ['15', 15.15],
        ['16', 33.382632],
        ['17', 40.22323],
        ['18', 19.11],
    '''

    for key, value in sorted(units.items()):
        if key in ['mm9', 'mm10']:
            coordinates = bp_coordinates
        else:
            coordinates = cm_coordinates

        for to_unit in sorted(value):
            req_dict = {
                'includeErrors': True,
                'includeInput': True,
                'includeExtra': True,
                'fromUnit': key,
                'toUnit': to_unit,
                'dataToConvert': coordinates
            }

            data = convert(req_dict)
            # iterate through all of the returned cordinates and print them out
            for res in data:
                print('{}\t{}\t{}'.format(key, to_unit, iterable_to_str(res, '\t')))

if __name__ == '__main__':
    main()