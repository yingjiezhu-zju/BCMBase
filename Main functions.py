from django.http import HttpResponse, FileResponse
from django.shortcuts import render, redirect, get_object_or_404
from subprocess import Popen, PIPE
import os
from django.conf import settings
from seawater_factor.models import SeawaterFactor
from soil_factor.models import SoilFactor
from math import radians, cos, sin, asin, sqrt
from django.db.models import Q
import csv
from datetime import datetime
from amplicon.models import Amplicon
from metagenome.models import Metagenome
from django.http import JsonResponse
from django.core.cache import cache
import logging
from microbe.models import Microbe
from genome.models import Genome
from builtins import map as builtin_map  
import subprocess 
import tempfile 
from django.urls import reverse 
from django import template

logger = logging.getLogger(__name__)

import re

def parse_coordinate(coord_str):
    if not coord_str:
        return None
    
    coord_str = coord_str.strip().upper()
    
    patterns = [
        r'^(-?\d+\.?\d*)\s*([NSEW]?)$', 
        r'^(\d+\.?\d*)\s*([NS])$',      
        r'^(\d+\.?\d*)\s*([EW])$',      
        r'^(-?\d+\.?\d*)$'             
    ]
    
    for pattern in patterns:
        match = re.match(pattern, coord_str)
        if match:
            value = float(match.group(1))
            direction = match.group(2) if len(match.groups()) > 1 else ''
            
           
            if direction in ['S', 'W']:
                value = -abs(value)
            elif direction in ['N', 'E']:
                value = abs(value)
            
            return value
    
    return None

def get_environment_factors(latitude, longitude, range_deg=0.5):

    if latitude is None or longitude is None:
        return {
            'seawater_factors': [],
            'soil_factors': [],
            'total_count': 0
        }
    

    lat_min = latitude - range_deg
    lat_max = latitude + range_deg
    lon_min = longitude - range_deg
    lon_max = longitude + range_deg
    

    seawater_factors = SeawaterFactor.objects.filter(
        latitude__gte=lat_min,
        latitude__lte=lat_max,
        longitude__gte=lon_min,
        longitude__lte=lon_max
    ).order_by('latitude', 'longitude')
    

    soil_factors = SoilFactor.objects.filter(
        latitude__gte=lat_min,
        latitude__lte=lat_max,
        longitude__gte=lon_min,
        longitude__lte=lon_max
    ).order_by('latitude', 'longitude')
    
    total_count = seawater_factors.count() + soil_factors.count()
    
    return {
        'seawater_factors': seawater_factors,
        'soil_factors': soil_factors,
        'total_count': total_count
    }

def parse_sample_coordinates(latitude_str, longitude_str):

    lat = parse_coordinate(latitude_str)
    lon = parse_coordinate(longitude_str)
    
    return lat, lon


BLAST_EXECUTABLE_PATH = '/usr/bin/blastn' 
BLAST_DB_PATH = '/home/zju/blast_db/my_nucleotide_db' 
TEMP_DIR = os.path.join(settings.BASE_DIR, 'temp_blast') 

os.makedirs(TEMP_DIR, exist_ok=True)

register = template.Library()

@register.filter
def split(value, delimiter):
    """Split a string by delimiter and return list."""
    return value.split(delimiter)

def home(request):
    return render(request, "home.html")

def network(request):
    return render(request, "network.html")

def haversine(lon1, lat1, lon2, lat2):

    try:
        coords = [float(lon1), float(lat1), float(lon2), float(lat2)]
      
        if not (-180 <= coords[0] <= 180 and -90 <= coords[1] <= 90 and
                -180 <= coords[2] <= 180 and -90 <= coords[3] <= 90):
            raise ValueError("Coordinates out of valid range")

        lon1, lat1, lon2, lat2 = list(builtin_map(radians, coords))

        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a))
        r = 6371  
        return c * r
    except (ValueError, TypeError) as e:
        logger.error(f"Error in haversine calculation: {str(e)}")
        return float('inf')  

def factors(request):

    context = {}

    search_type = request.GET.get('search_type')

    longitude = request.GET.get('longitude')
    latitude = request.GET.get('latitude')
    distance = float(request.GET.get('distance', 500))
    factor_type = request.GET.get('type')

    context.update({
        'distance': distance,
        'factor_type': factor_type,
        'search_type': search_type
    })

    if search_type:
        try:
            if search_type == 'seawater':

                salinity_min = request.GET.get('salinity_min')
                salinity_max = request.GET.get('salinity_max')
                temp_min = request.GET.get('temp_min')
                temp_max = request.GET.get('temp_max')
                ph_min = request.GET.get('ph_min')
                ph_max = request.GET.get('ph_max')

                query = Q()
                if salinity_min and salinity_max:
                    query &= Q(salinity__gte=float(salinity_min)) & Q(salinity__lte=float(salinity_max))
                if temp_min and temp_max:
                    query &= Q(ocean_temperature__gte=float(temp_min)) & Q(ocean_temperature__lte=float(temp_max))
                if ph_min and ph_max:
                    query &= Q(ph__gte=float(ph_min)) & Q(ph__lte=float(ph_max))

                seawater_factors = list(SeawaterFactor.objects.filter(query))

                logger.debug(f"Found {len(seawater_factors)} seawater factors")
                for factor in seawater_factors:
                    logger.debug(f"Seawater factor at ({factor.latitude}, {factor.longitude})")

                context['seawater_factors'] = seawater_factors
                context['searched'] = True

            elif search_type == 'soil':
                organic_carbon_min = request.GET.get('organic_carbon_min')
                organic_carbon_max = request.GET.get('organic_carbon_max')
                total_carbon_min = request.GET.get('total_carbon_min')
                total_carbon_max = request.GET.get('total_carbon_max')
                organic_matter_min = request.GET.get('organic_matter_min')
                organic_matter_max = request.GET.get('organic_matter_max')

                query = Q()
                if organic_carbon_min and organic_carbon_max:
                    query &= Q(average_organic_carbon__gte=float(organic_carbon_min)) & Q(average_organic_carbon__lte=float(organic_carbon_max))
                if total_carbon_min and total_carbon_max:
                    query &= Q(average_total_carbon__gte=float(total_carbon_min)) & Q(average_total_carbon__lte=float(total_carbon_max))
                if organic_matter_min and organic_matter_max:
                    query &= Q(average_organic_matter__gte=float(organic_matter_min)) & Q(average_organic_matter__lte=float(organic_matter_max))

                soil_factors = list(SoilFactor.objects.filter(query))

                logger.debug(f"Found {len(soil_factors)} soil factors")
                for factor in soil_factors:
                    logger.debug(f"Soil factor at ({factor.latitude}, {factor.longitude})")

                context['soil_factors'] = soil_factors
                context['searched'] = True

        except ValueError as e:
            logger.error(f"Parameter error: {str(e)}")
            context['error'] = "Invalid parameter values. Please check your input."

    elif longitude and latitude:
        try:
            longitude = float(longitude)
            latitude = float(latitude)

            if not (-180 <= longitude <= 180 and -90 <= latitude <= 90):
                raise ValueError("Coordinates out of range")

            context['searched'] = True

            if factor_type != 'soil':
                seawater_factors = list(SeawaterFactor.objects.all())
                for factor in seawater_factors:
                    factor.distance = haversine(longitude, latitude, factor.longitude, factor.latitude)

                seawater_factors = sorted(
                    [f for f in seawater_factors if f.distance <= distance],
                    key=lambda x: x.distance
                )
                context['seawater_factors'] = seawater_factors

            if factor_type != 'seawater':
                soil_factors = list(SoilFactor.objects.all())
                for factor in soil_factors:
                    factor.distance = haversine(longitude, latitude, factor.longitude, factor.latitude)
                soil_factors = sorted(
                    [f for f in soil_factors if f.distance <= distance],
                    key=lambda x: x.distance
                )
                context['soil_factors'] = soil_factors

        except ValueError as e:
            context['error'] = "Invalid coordinates. Please ensure longitude is between -180 and 180, and latitude is between -90 and 90."
            print(f"Coordinate error in factors view: {str(e)}")

    return render(request, "factors.html", context)

def seawater_detail(request, factor_id):

    factor = get_object_or_404(SeawaterFactor, id=factor_id)
    return render(request, "seawater_detail.html", {'factor': factor})

def soil_detail(request, factor_id):

    factor = get_object_or_404(SoilFactor, id=factor_id)
    return render(request, "soil_detail.html", {'factor': factor})

def export_seawater_factors(request):

    try:
        longitude = float(request.GET.get('longitude'))
        latitude = float(request.GET.get('latitude'))
        distance = float(request.GET.get('distance', 500))

        seawater_factors = list(SeawaterFactor.objects.all())
        filtered_factors = []
        for factor in seawater_factors:
            factor.distance = haversine(longitude, latitude, factor.longitude, factor.latitude)
            if factor.distance <= distance:
                filtered_factors.append(factor)
        filtered_factors.sort(key=lambda x: x.distance)
        response = HttpResponse(content_type='text/tab-separated-values')
        response['Content-Disposition'] = f'attachment; filename="seawater_factors_{datetime.now().strftime("%Y%m%d_%H%M%S")}.tsv"'
        writer = csv.writer(response, delimiter='\t')
        writer.writerow([
            'ID', 'Latitude', 'Longitude', 'Distance (km)',
            'Salinity (PSU)',
            'Silicate (μmol/L)',
            'Phosphate (μmol/L)',
            'Nitrate (μmol/L)',
            'Iron (μmol/L)',
            'pH',
            'Dissolved Oxygen (μmol/kg)',
            'Temperature (°C)',
            'Seawater Direction (degrees)',
            'Seawater Speed (m/s)',
            'Primary Productivity (mg C/m³/day)',
            'Created At', 'Updated At'
        ])

        for factor in filtered_factors:
            writer.writerow([
                factor.id,
                factor.latitude,
                factor.longitude,
                f"{factor.distance:.2f}",
                factor.salinity,
                factor.silicate,
                factor.phosphate,
                factor.nitrate,
                factor.iron,
                factor.ph,
                factor.dissolved_oxygen,
                factor.ocean_temperature,
                factor.seawater_direction,
                factor.seawater_speed,
                factor.primary_productivity,
                factor.created_at.strftime("%Y-%m-%d %H:%M:%S"),
                factor.updated_at.strftime("%Y-%m-%d %H:%M:%S")
            ])

        return response

    except (ValueError, TypeError) as e:
        return HttpResponse(f"Error: {str(e)}", status=400)

def export_soil_factors(request):
    """
    Export soil factors search results to TSV based on the SoilFactor model fields
    """
    try:
        longitude = float(request.GET.get('longitude'))
        latitude = float(request.GET.get('latitude'))
        distance = float(request.GET.get('distance', 500))
        soil_factors = list(SoilFactor.objects.all())
        filtered_factors = []
        for factor in soil_factors:
            factor.distance = haversine(longitude, latitude, factor.longitude, factor.latitude)
            if factor.distance <= distance:
                filtered_factors.append(factor)

        filtered_factors.sort(key=lambda x: x.distance)
        response = HttpResponse(content_type='text/tab-separated-values')
        response['Content-Disposition'] = f'attachment; filename="soil_factors_{datetime.now().strftime("%Y%m%d_%H%M%S")}.tsv"'
        writer = csv.writer(response, delimiter='\t')
        writer.writerow([
            'ID', 'Latitude', 'Longitude', 'Distance (km)',
            'Average Organic Carbon (%)',
            'Organic Carbon Upper Depth (cm)', 'Organic Carbon Lower Depth (cm)',
            'Average Total Carbon (%)',
            'Total Carbon Upper Depth (cm)', 'Total Carbon Lower Depth (cm)',
            'Average Organic Matter (%)',
            'Organic Matter Upper Depth (cm)', 'Organic Matter Lower Depth (cm)',
            'Created At', 'Updated At'
        ])

        for factor in filtered_factors:
            writer.writerow([
                factor.id,
                factor.latitude,
                factor.longitude,
                f"{factor.distance:.2f}",
                factor.average_organic_carbon,
                factor.organic_carbon_upper_depth,
                factor.organic_carbon_lower_depth,
                factor.average_total_carbon,
                factor.total_carbon_upper_depth,
                factor.total_carbon_lower_depth,
                factor.average_organic_matter,
                factor.organic_matter_upper_depth,
                factor.organic_matter_lower_depth,
                factor.created_at.strftime("%Y-%m-%d %H:%M:%S"),
                factor.updated_at.strftime("%Y-%m-%d %H:%M:%S")
            ])

        return response

    except (ValueError, TypeError) as e:
        return HttpResponse(f"Error: {str(e)}", status=400)

def map_view(request):
    return render(request, "tools/map.html")

def map_data(request):
    try:
        logger.info(f"Received map data request with parameters: {request.GET}")

        sample_types = request.GET.getlist('sample_types', ['amplicon', 'metagenome'])
        date_from = request.GET.get('date_from', '').strip()
        date_to = request.GET.get('date_to', '').strip()
        host = request.GET.get('host', '').strip().lower()  
        habitat = request.GET.get('habitat', '').strip().lower()
        location = request.GET.get('location', '').strip().lower()

        valid_types = ['amplicon', 'metagenome']
        sample_types = [t for t in sample_types if t in valid_types]
        if not sample_types:
            return JsonResponse({
                'status': 'error',
                'message': 'Invalid sample types'
            }, status=400)

        logger.info(f"""Filter parameters:
            Sample types: {sample_types}
            Date range: {date_from} - {date_to}
            Host: {host}
            Habitat: {habitat}
            Location: {location}
        """)

        if date_from or date_to:
            logger.info("Validating date format")
            try:
                if date_from:
                    datetime.strptime(date_from, '%Y-%m-%d')
                if date_to:
                    datetime.strptime(date_to, '%Y-%m-%d')
                    if date_from and date_to < date_from:
                        return JsonResponse({
                            'status': 'error',
                            'message': 'End date cannot be earlier than start date'
                        }, status=400)
            except ValueError as e:
                logger.error(f"Date validation error: {str(e)}")
                return JsonResponse({
                    'status': 'error',
                    'message': 'Invalid date format. Use YYYY-MM-DD'
                }, status=400)

        microbe_rank = request.GET.get('microbe_rank', '').strip()
        microbe_term = request.GET.get('microbe_term', '').strip().lower()
        microbe_function = request.GET.get('microbe_function', '').strip().lower()
        sample_runs = set()
        if microbe_rank and microbe_term:
            microbe_query = Q()

            if microbe_rank == 'species':
                microbe_query &= Q(species__icontains=microbe_term)
            elif microbe_rank == 'genus':
                microbe_query &= Q(genus__icontains=microbe_term)
            elif microbe_rank == 'family':
                microbe_query &= Q(family__icontains=microbe_term)

            if microbe_function:
                microbe_query &= Q(functions__contains=[microbe_function])

            microbes = Microbe.objects.filter(microbe_query)
            for microbe in microbes:
                sample_runs.update(microbe.sources)
        data = []
        if 'amplicon' in sample_types:
            logger.info("Processing Amplicon data")
            try:
                amplicon_query = Q()

                if date_from:
                    amplicon_query &= Q(collection_date__gte=date_from)
                if date_to:
                    amplicon_query &= Q(collection_date__lte=date_to)
                if host:
                    amplicon_query &= Q(host__icontains=host)
                if habitat:
                    amplicon_query &= Q(habitat__icontains=habitat)
                if location:
                    amplicon_query &= Q(geo_location__icontains=location)

                if sample_runs:
                    amplicon_query &= Q(run__in=sample_runs)

                logger.info(f"Amplicon query: {amplicon_query}")

                amplicon_samples = (Amplicon.objects.filter(amplicon_query)
                                  .exclude(latitude='')
                                  .exclude(longitude='')
                                  .only('latitude', 'longitude', 'run', 'host',
                                      'geo_location', 'collection_date', 'habitat'))

                logger.info(f"Found {amplicon_samples.count()} Amplicon samples")

                for sample in amplicon_samples:
                    try:
                        lat = float(sample.latitude.replace('N', '').replace('S', '').strip())
                        lng = float(sample.longitude.replace('E', '').replace('W', '').strip())

                        if 'S' in sample.latitude:
                            lat = -lat
                        if 'W' in sample.longitude:
                            lng = -lng

                        if not (-90 <= lat <= 90 and -180 <= lng <= 180):
                            logger.warning(f"Invalid coordinates for sample {sample.run}: {lat}, {lng}")
                            continue

                        data.append({
                            'type': 'amplicon',
                            'lat': lat,
                            'lng': lng,
                            'run': sample.run,
                            'host': sample.host,
                            'country': sample.geo_location.split(':')[0] if sample.geo_location else None,
                            'collection_date': sample.collection_date,
                            'habitat': sample.habitat,
                            'url': f'/amplicon/{sample.run}/'
                        })
                    except (ValueError, AttributeError) as e:
                        logger.warning(f"Error processing amplicon sample {sample.run}: {str(e)}")
                        continue

            except Exception as e:
                logger.error(f"Error processing Amplicon data: {str(e)}")

        if 'metagenome' in sample_types:
            logger.info("Processing Metagenome data")
            try:
                metagenome_query = Q()
                if date_from:
                    metagenome_query &= Q(collection_date__gte=date_from)
                if date_to:
                    metagenome_query &= Q(collection_date__lte=date_to)
                if host:
                    metagenome_query &= Q(host__icontains=host)
                if habitat:
                    metagenome_query &= Q(habitat__icontains=habitat)
                if location:
                    metagenome_query &= Q(geo_location__icontains=location)

                if sample_runs:
                    metagenome_query &= Q(run_id__in=sample_runs)

                logger.info(f"Metagenome query: {metagenome_query}")

                metagenome_samples = (Metagenome.objects.filter(metagenome_query)
                                    .exclude(latitude='')
                                    .exclude(longitude='')
                                    .only('latitude', 'longitude', 'run_id', 'host',
                                        'geo_location', 'collection_date', 'habitat'))

                logger.info(f"Found {metagenome_samples.count()} Metagenome samples")

                for sample in metagenome_samples:
                    try:
                        lat = float(sample.latitude.replace('N', '').replace('S', '').strip())
                        lng = float(sample.longitude.replace('E', '').replace('W', '').strip())

                        if 'S' in sample.latitude:
                            lat = -lat
                        if 'W' in sample.longitude:
                            lng = -lng

                        if not (-90 <= lat <= 90 and -180 <= lng <= 180):
                            logger.warning(f"Invalid coordinates for sample {sample.run_id}: {lat}, {lng}")
                            continue

                        data.append({
                            'type': 'metagenome',
                            'lat': lat,
                            'lng': lng,
                            'run': sample.run_id,
                            'host': sample.host,
                            'country': sample.geo_location.split(':')[0] if sample.geo_location else None,
                            'collection_date': sample.collection_date,
                            'habitat': sample.habitat,
                            'url': f'/metagenome/{sample.run_id}/'
                        })
                    except (ValueError, AttributeError) as e:
                        logger.warning(f"Error processing metagenome sample {sample.run_id}: {str(e)}")
                        continue

            except Exception as e:
                logger.error(f"Error processing Metagenome data: {str(e)}")
                # 继续处理其他数据类型，不中断整个请求

        # 准备响应数据
        response_data = {
            'status': 'success',
            'data': data,
            'total_count': len(data),
            'amplicon_count': sum(1 for x in data if x['type'] == 'amplicon'),
            'metagenome_count': sum(1 for x in data if x['type'] == 'metagenome')
        }

        logger.info(f"""Response summary:
            Total samples: {response_data['total_count']}
            Amplicon samples: {response_data['amplicon_count']}
            Metagenome samples: {response_data['metagenome_count']}
        """)

        return JsonResponse(response_data)

    except Exception as e:
        logger.error(f"Error in map_data view: {str(e)}", exc_info=True)
        return JsonResponse({
            'status': 'error',
            'message': str(e)
        }, status=500)

def factors_parameters(request):
    """
    View function for searching environmental factors by parameters.
    """
    context = {
        'search_type': request.GET.get('search_type', 'seawater')
    }

    if request.GET:
        try:
            if context['search_type'] == 'seawater':
                silicate_min = request.GET.get('silicate_min')
                silicate_max = request.GET.get('silicate_max')
                phosphate_min = request.GET.get('phosphate_min')
                phosphate_max = request.GET.get('phosphate_max')
                nitrate_min = request.GET.get('nitrate_min')
                nitrate_max = request.GET.get('nitrate_max')
                iron_min = request.GET.get('iron_min')
                iron_max = request.GET.get('iron_max')
                dissolved_oxygen_min = request.GET.get('dissolved_oxygen_min')
                dissolved_oxygen_max = request.GET.get('dissolved_oxygen_max')
                primary_productivity_min = request.GET.get('primary_productivity_min')
                primary_productivity_max = request.GET.get('primary_productivity_max')
                query = Q()
                if silicate_min and silicate_max:
                    query &= Q(silicate__gte=float(silicate_min)) & Q(silicate__lte=float(silicate_max))
                if phosphate_min and phosphate_max:
                    query &= Q(phosphate__gte=float(phosphate_min)) & Q(phosphate__lte=float(phosphate_max))
                if nitrate_min and nitrate_max:
                    query &= Q(nitrate__gte=float(nitrate_min)) & Q(nitrate__lte=float(nitrate_max))
                if iron_min and iron_max:
                    query &= Q(iron__gte=float(iron_min)) & Q(iron__lte=float(iron_max))
                if dissolved_oxygen_min and dissolved_oxygen_max:
                    query &= Q(dissolved_oxygen__gte=float(dissolved_oxygen_min)) & Q(dissolved_oxygen__lte=float(dissolved_oxygen_max))
                if primary_productivity_min and primary_productivity_max:
                    query &= Q(primary_productivity__gte=float(primary_productivity_min)) & Q(primary_productivity__lte=float(primary_productivity_max))
                salinity_min = request.GET.get('salinity_min')
                salinity_max = request.GET.get('salinity_max')
                temp_min = request.GET.get('temp_min')
                temp_max = request.GET.get('temp_max')
                ph_min = request.GET.get('ph_min')
                ph_max = request.GET.get('ph_max')
                if salinity_min and salinity_max:
                    query &= Q(salinity__gte=float(salinity_min)) & Q(salinity__lte=float(salinity_max))
                if temp_min and temp_max:
                    query &= Q(ocean_temperature__gte=float(temp_min)) & Q(ocean_temperature__lte=float(temp_max))
                if ph_min and ph_max:
                    query &= Q(ph__gte=float(ph_min)) & Q(ph__lte=float(ph_max))
                logger.debug(f"Seawater query: {query}")
                seawater_factors = SeawaterFactor.objects.filter(query)
                total_count = seawater_factors.count()
                logger.debug(f"Found {total_count} seawater factors")

                if total_count > 5000:
                    context.update({
                        'searched': True,
                        'too_many_results': True,
                        'total_count': total_count
                    })
                else:
                    factors_list = list(seawater_factors)
                    for factor in factors_list[:5]:  
                        logger.debug(f"Factor coordinates: lat={factor.latitude}, lon={factor.longitude}")

                    context.update({
                        'searched': True,
                        'seawater_factors': factors_list
                    })

            elif context['search_type'] == 'soil':
                organic_carbon_min = request.GET.get('organic_carbon_min')
                organic_carbon_max = request.GET.get('organic_carbon_max')
                organic_carbon_upper_depth = request.GET.get('organic_carbon_upper_depth')
                organic_carbon_lower_depth = request.GET.get('organic_carbon_lower_depth')

                total_carbon_min = request.GET.get('total_carbon_min')
                total_carbon_max = request.GET.get('total_carbon_max')
                total_carbon_upper_depth = request.GET.get('total_carbon_upper_depth')
                total_carbon_lower_depth = request.GET.get('total_carbon_lower_depth')

                organic_matter_min = request.GET.get('organic_matter_min')
                organic_matter_max = request.GET.get('organic_matter_max')
                organic_matter_upper_depth = request.GET.get('organic_matter_upper_depth')
                organic_matter_lower_depth = request.GET.get('organic_matter_lower_depth')

                query = Q()
                if organic_carbon_min and organic_carbon_max:
                    query &= Q(average_organic_carbon__gte=float(organic_carbon_min)) & Q(average_organic_carbon__lte=float(organic_carbon_max))
                if organic_carbon_upper_depth:
                    query &= Q(organic_carbon_upper_depth=int(organic_carbon_upper_depth))
                if organic_carbon_lower_depth:
                    query &= Q(organic_carbon_lower_depth=int(organic_carbon_lower_depth))
                if total_carbon_min and total_carbon_max:
                    query &= Q(average_total_carbon__gte=float(total_carbon_min)) & Q(average_total_carbon__lte=float(total_carbon_max))
                if total_carbon_upper_depth:
                    query &= Q(total_carbon_upper_depth=int(total_carbon_upper_depth))
                if total_carbon_lower_depth:
                    query &= Q(total_carbon_lower_depth=int(total_carbon_lower_depth))
                if organic_matter_min and organic_matter_max:
                    query &= Q(average_organic_matter__gte=float(organic_matter_min)) & Q(average_organic_matter__lte=float(organic_matter_max))
                if organic_matter_upper_depth:
                    query &= Q(organic_matter_upper_depth=int(organic_matter_upper_depth))
                if organic_matter_lower_depth:
                    query &= Q(organic_matter_lower_depth=int(organic_matter_lower_depth))
                logger.debug(f"Soil query: {query}")
                soil_factors = SoilFactor.objects.filter(query)
                total_count = soil_factors.count()
                logger.debug(f"Found {total_count} soil factors")

                if total_count > 5000:
                    context.update({
                        'searched': True,
                        'too_many_results': True,
                        'total_count': total_count
                    })
                else:
                    factors_list = list(soil_factors)
                    for factor in factors_list[:5]:  
                        logger.debug(f"Factor coordinates: lat={factor.latitude}, lon={factor.longitude}")

                    context.update({
                        'searched': True,
                        'soil_factors': factors_list
                    })

        except ValueError as e:
            context['error'] = "Invalid parameter values. Please check your input."
            logger.error(f"Parameter error: {str(e)}")
    logger.debug(f"Context keys: {context.keys()}")
    if 'seawater_factors' in context:
        logger.debug(f"Seawater factors count in context: {len(context['seawater_factors'])}")
    if 'soil_factors' in context:
        logger.debug(f"Soil factors count in context: {len(context['soil_factors'])}")

    return render(request, "factors_parameters.html", context)

def blast_tool(request):
    """
    View function for the Diamond BLAST tool page.
    Handles displaying the form (GET) and processing submission (POST).
    """
    DIAMOND_EXECUTABLE_PATH = '/usr/local/bin/diamond'  

    AVAILABLE_DATABASES = {
        'mangrove': {
            'path': '/home/zju/GBCBase/media/Mangrove.gene.dmnd',
            'name': 'Mangrove Gene Database',
            'description': 'Microbial gene catalog for mangrove ecosystems, including gene sequences unique to mangrove environments.'
        },
        'saltmarsh': {
            'path': '/home/zju/GBCBase/media/Saltmarsh.gene.dmnd',
            'name': 'Saltmarsh Gene Database',
            'description': 'Microbial gene catalog for salt marsh wetland ecosystems, suitable for research on halotolerant microbes.'
        },
        'seagrass': {
            'path': '/home/zju/GBCBase/media/Seagrass.gene.dmnd',
            'name': 'Seagrass Gene Database',
            'description': 'Microbial gene catalog for seagrass bed ecosystems, covering microbes from the seagrass rhizosphere and sediments.'
        }
    }

    if request.method == 'POST':
        sequence_input = request.POST.get('sequence', '').strip()
        sequence_file = request.FILES.get('sequence_file')
        blast_program = request.POST.get('blast_program', 'blastp')
        evalue = request.POST.get('evalue', '1e-5')
        max_target_seqs = request.POST.get('max_target_seqs', '5')
        database_choice = request.POST.get('database', 'mangrove') 

        if database_choice not in AVAILABLE_DATABASES:
            database_choice = 'mangrove' 

        selected_db = AVAILABLE_DATABASES[database_choice]
        DIAMOND_DB_PATH = selected_db['path']

        logger.info(f"Using Diamond database: {selected_db['name']} at {DIAMOND_DB_PATH}")

        if not sequence_input and not sequence_file:
            return render(request, 'tools/blast.html', {
                'error': 'Please provide a sequence or upload a file.',
                'blast_program': blast_program,
                'evalue': evalue,
                'max_target_seqs': max_target_seqs,
                'database': database_choice,
                'available_databases': AVAILABLE_DATABASES
            })

        try:

            os.makedirs(os.path.dirname(TEMP_DIR), exist_ok=True)

            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa", dir=TEMP_DIR) as query_file:

                if sequence_file:
                    content = sequence_file.read().decode('utf-8')
                    query_file.write(content)

                else:
                    if not sequence_input.startswith('>'):
                        query_file.write(">Query\n")
                    query_file.write(sequence_input)

                query_filepath = query_file.name

            output_file = tempfile.NamedTemporaryFile(delete=False, suffix=".out", dir=TEMP_DIR)
            output_filepath = output_file.name
            output_file.close()

            diamond_command = [
                DIAMOND_EXECUTABLE_PATH,
                blast_program, 
                '--query', query_filepath,
                '--db', DIAMOND_DB_PATH,
                '--out', output_filepath,
                '--max-target-seqs', max_target_seqs,
                '--evalue', evalue,
                '--outfmt', '6',  
                '--threads', '2'  
            ]

            logger.info(f"Running Diamond command: {' '.join(diamond_command)}")

            process = subprocess.Popen(
                diamond_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            stdout, stderr = process.communicate(timeout=120)  

            if process.returncode != 0:
                logger.error(f"Diamond error: {stderr}")
                return render(request, 'tools/blast.html', {
                    'error': f"Diamond BLAST search failed: {stderr[:200]}...",
                    'sequence': sequence_input,
                    'blast_program': blast_program,
                    'evalue': evalue,
                    'max_target_seqs': max_target_seqs,
                    'database': database_choice,
                    'available_databases': AVAILABLE_DATABASES
                })

            with open(output_filepath, 'r') as f:
                blast_results = f.readlines()
            enriched_results = []
            column_names = [
                'query_id', 'subject_id', 'percent_identity', 'alignment_length',
                'mismatches', 'gap_opens', 'query_start', 'query_end',
                'subject_start', 'subject_end', 'evalue', 'bit_score'
            ]

            for line in blast_results:
                values = line.strip().split('\t')
                if len(values) == len(column_names):
                    result_dict = dict(zip(column_names, values))
                    query_parts = result_dict['query_id'].split('|')
                    query_info = {
                        'accession': query_parts[0] if len(query_parts) > 0 else "",
                        'type': query_parts[1] if len(query_parts) > 1 else "",
                        'description': query_parts[2] if len(query_parts) > 2 else ""
                    }

                    subject_parts = result_dict['subject_id'].split('|')
                    subject_info = {
                        'source': subject_parts[0] if len(subject_parts) > 0 else "",
                        'cog_category': subject_parts[1] if len(subject_parts) > 1 else "",
                        'description': (subject_parts[2] if len(subject_parts) > 2 else "").replace("_", " "),
                        'preferred_name': (subject_parts[3] if len(subject_parts) > 3 else "").replace("_", " "),
                        'cazy': subject_parts[4] if len(subject_parts) > 4 else "-"
                    }
                    identity_value = float(result_dict['percent_identity'])
                    if identity_value >= 90:
                        identity_class = "bg-green-500"
                    elif identity_value >= 70:
                        identity_class = "bg-blue-500"
                    elif identity_value >= 50:
                        identity_class = "bg-yellow-500"
                    else:
                        identity_class = "bg-red-500"

                    evalue = result_dict['evalue']
                    if evalue == '0.0':
                        evalue_class = "text-green-600 font-semibold"
                    elif 'e-' in evalue:
                        try:
                            e_power = int(evalue.split('e-')[1])
                            if e_power > 100:
                                evalue_class = "text-green-600 font-semibold"
                            elif e_power > 50:
                                evalue_class = "text-blue-600"
                            elif e_power > 10:
                                evalue_class = "text-gray-700"
                            else:
                                evalue_class = "text-gray-500"
                        except (IndexError, ValueError):
                            evalue_class = "text-gray-500"
                    else:
                        evalue_class = "text-gray-500"

                    bit_score = float(result_dict['bit_score'])
                    if bit_score > 1000:
                        bit_score_class = "text-green-600 font-semibold"
                    elif bit_score > 500:
                        bit_score_class = "text-blue-600"
                    elif bit_score > 200:
                        bit_score_class = "text-gray-700"
                    else:
                        bit_score_class = "text-gray-500"
                    subject_type_color = "bg-blue-100 text-blue-800"  
                    result_dict.update({
                        'query_info': query_info,
                        'subject_info': subject_info,
                        'identity_class': identity_class,
                        'evalue_class': evalue_class,
                        'bit_score_class': bit_score_class,
                        'subject_type_color': subject_type_color
                    })

                    enriched_results.append(result_dict)

            request.session['blast_raw_result'] = ''.join(blast_results)
            original_result_count = len(enriched_results)
            filtered_results = [r for r in enriched_results if float(r['percent_identity']) >= 50]
            filtered_count = original_result_count - len(filtered_results)

            request.session['blast_parameters'] = {
                'program': blast_program,
                'evalue': evalue,
                'max_target_seqs': max_target_seqs,
                'database': selected_db['name'],
                'database_description': selected_db['description'],
                'original_count': original_result_count,
                'filtered_count': filtered_count
            }

            os.remove(query_filepath)
            os.remove(output_filepath)
            return render(request, 'tools/blast_result.html', {
                'results': filtered_results,
                'parameters': request.session['blast_parameters'],
                'result_count': len(filtered_results),
                'has_filtered': filtered_count > 0,
                'filtered_count': filtered_count
            })

        except subprocess.TimeoutExpired:
            logger.error("Diamond command timed out")

            if 'query_filepath' in locals() and os.path.exists(query_filepath):
                os.remove(query_filepath)
            if 'output_filepath' in locals() and os.path.exists(output_filepath):
                os.remove(output_filepath)

            request.session['timeout_sequence'] = sequence_input
            request.session['timeout_parameters'] = {
                'blast_program': blast_program,
                'evalue': evalue,
                'max_target_seqs': max_target_seqs,
                'database': database_choice,
                'database_name': selected_db['name']
            }

            return redirect('blast_timeout')

        except Exception as e:
            logger.error(f"An error occurred during Diamond execution: {str(e)}", exc_info=True)

            if 'query_filepath' in locals() and os.path.exists(query_filepath):
                os.remove(query_filepath)
            if 'output_filepath' in locals() and os.path.exists(output_filepath):
                os.remove(output_filepath)

            return render(request, 'tools/blast.html', {
                'error': f'An unexpected error occurred: {str(e)}',
                'sequence': sequence_input,
                'blast_program': blast_program,
                'evalue': evalue,
                'max_target_seqs': max_target_seqs,
                'database': database_choice,
                'available_databases': AVAILABLE_DATABASES
            })

    return render(request, 'tools/blast.html', {
        'available_databases': AVAILABLE_DATABASES
    })

def blast_timeout(request):
    timeout_sequence = request.session.get('timeout_sequence', '')
    timeout_parameters = request.session.get('timeout_parameters', {})
    sequence_length = 0
    if timeout_sequence:
        lines = timeout_sequence.strip().split('\n')
        sequence_data = ''.join([line for line in lines if not line.startswith('>')])
        sequence_length = len(sequence_data.replace(' ', '').replace('\t', ''))

    context = {
        'sequence_length': sequence_length,
        'parameters': timeout_parameters,
        'has_sequence_data': bool(timeout_sequence)
    }

    if 'timeout_sequence' in request.session:
        del request.session['timeout_sequence']
    if 'timeout_parameters' in request.session:
        del request.session['timeout_parameters']

    return render(request, 'tools/blast_timeout.html', context)

def export_blast_result(request):
    """
    Export the BLAST results as a TSV file
    """
    raw_result = request.session.get('blast_raw_result')
    if not raw_result:
        return HttpResponse("No BLAST result found to export.", status=404)

    response = HttpResponse(raw_result, content_type='text/tab-separated-values')
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    response['Content-Disposition'] = f'attachment; filename="diamond_blast_result_{timestamp}.tsv"'

    return response

def help_batch(request):
    """
    View function for the batch search help page
    """
    return render(request, 'help_batch.html')

def help_interactive_map(request):
    return render(request, 'help_interactive_map.html')

def help_catalog(request):
    context = {
        'section': 'catalog'  
    }
    return render(request, 'help_catalog.html', context)

def help_network(request):
    context = {
        'section': 'network'  
    }
    return render(request, 'help_network.html', context)

def help_bcvhp(request):
    """
    View function for the BCVHP help page
    """
    context = {
        'section': 'bcvhp'
    }
    return render(request, 'help_bcvhp.html', context)

def handler404(request, exception):
    """Custom 404 error handler"""
    return render(request, 'errors/404.html', status=404)

def handler500(request):
    """Custom 500 error handler"""
    return render(request, 'errors/500.html', status=500)

def handler504(request, exception=None):

    if 'blast' in request.path:
        return redirect('blast_timeout')
    return render(request, 'errors/504.html', status=504)

def bcvhp_view(request):

    available_databases = {}  
    return render(request, 'tools/bcvhp.html', {
        'available_databases': available_databases
    })

BCVHP_BLASTN_PATH = '/software/ncbi-blast-2.6.0+/bin/blastn'
BCVHP_CRISPR_DB = '/home/zju/GBCBase/media/bcvhp/bcvhp_crispr/bcvhp_crispr.blastdb'
BCVHP_TRNA_DB = '/home/zju/GBCBase/media/bcvhp/bcvhp_trna/bcvhp_trna'
BCVHP_TEMP_DIR = os.path.join(settings.BASE_DIR, 'temp_bcvhp')
os.makedirs(BCVHP_TEMP_DIR, exist_ok=True)

def bcvhp_tool(request):

    available_databases = {
        'bcvhp_crispr': {
            'path': BCVHP_CRISPR_DB,
            'name': 'BCVHP CRISPR Database',
            'description': 'CRISPR spacer database for virus host prediction.'
        },
        'bcvhp_trna': {
            'path': BCVHP_TRNA_DB,
            'name': 'BCVHP tRNA Database',
            'description': 'tRNA database for virus host prediction based on tRNA sequences.'
        }

    }
    
    if request.method == 'POST':
        sequence_input = request.POST.get('sequence', '').strip()
        sequence_file = request.FILES.get('sequence_file')
        blast_program = request.POST.get('blast_program', 'blastn')
        evalue = request.POST.get('evalue', '1e-5')
        max_target_seqs = request.POST.get('max_target_seqs', '5')
        culling_limit = request.POST.get('culling_limit', '1')
        database_choice = request.POST.get('database', 'bcvhp_crispr')

        if database_choice not in available_databases:
            database_choice = 'bcvhp_crispr'
        selected_db = available_databases[database_choice]
        db_path = selected_db['path']

        if not sequence_input and not sequence_file:
            return render(request, 'tools/bcvhp.html', {
                'error': 'Please provide a sequence or upload a file.',
                'blast_program': blast_program,
                'evalue': evalue,
                'max_target_seqs': max_target_seqs,
                'culling_limit': culling_limit,
                'database': database_choice,
                'available_databases': available_databases
            })

        try:
            with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa", dir=BCVHP_TEMP_DIR) as query_file:
                if sequence_file:
                    content = sequence_file.read().decode('utf-8')
                    query_file.write(content)
                else:
                    if not sequence_input.startswith('>'):
                        query_file.write(">Query\n")
                    query_file.write(sequence_input)
                query_filepath = query_file.name

            output_file = tempfile.NamedTemporaryFile(delete=False, suffix=".out", dir=BCVHP_TEMP_DIR)
            output_filepath = output_file.name
            output_file.close()

            blastn_command = [
                BCVHP_BLASTN_PATH,
                '-query', query_filepath,
                '-db', db_path,
                '-out', output_filepath,
                '-outfmt', '6',
                '-max_target_seqs', max_target_seqs,
                '-evalue', evalue,
                '-num_threads', '1',
                '-perc_identity', '95',
                '-culling_limit', culling_limit
            ]

            process = subprocess.Popen(
                blastn_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            stdout, stderr = process.communicate(timeout=120)

            if process.returncode != 0:
                return render(request, 'tools/bcvhp.html', {
                    'error': f"BLASTN search failed: {stderr[:200]}...",
                    'sequence': sequence_input,
                    'blast_program': blast_program,
                    'evalue': evalue,
                    'max_target_seqs': max_target_seqs,
                    'culling_limit': culling_limit,
                    'database': database_choice,
                    'available_databases': available_databases
                })

            with open(output_filepath, 'r') as f:
                blast_results = f.readlines()
            results = []
            for line in blast_results:
                values = line.strip().split('\t')
                if len(values) >= 12:
                    identity_value = float(values[2])
                    if identity_value >= 90:
                        identity_class = "bg-green-500"
                    elif identity_value >= 70:
                        identity_class = "bg-blue-500"
                    elif identity_value >= 50:
                        identity_class = "bg-yellow-500"
                    else:
                        identity_class = "bg-red-500"
                    evalue = values[10]
                    if evalue == '0.0':
                        evalue_class = "text-green-600 font-semibold"
                    elif 'e-' in evalue:
                        try:
                            e_power = int(evalue.split('e-')[1])
                            if e_power > 100:
                                evalue_class = "text-green-600 font-semibold"
                            elif e_power > 50:
                                evalue_class = "text-blue-600"
                            elif e_power > 10:
                                evalue_class = "text-gray-700"
                            else:
                                evalue_class = "text-gray-500"
                        except (IndexError, ValueError):
                            evalue_class = "text-gray-500"
                    else:
                        evalue_class = "text-gray-500"

                    bit_score = float(values[11])
                    if bit_score > 1000:
                        bit_score_class = "text-green-600 font-semibold"
                    elif bit_score > 500:
                        bit_score_class = "text-blue-600"
                    elif bit_score > 200:
                        bit_score_class = "text-gray-700"
                    else:
                        bit_score_class = "text-gray-500"

                    host_genome_id = values[1].split('|')[0] if '|' in values[1] else values[1]
                    annotation = '-'
                    gtdb_detail = ''
                    try:
                        genome = Genome.objects.get(bin_id=host_genome_id)
                        if genome.marker_lineage:
                            annotation = genome.marker_lineage
                        if genome.gtdb_classification:
                            gtdb_detail = genome.gtdb_classification
                    except Genome.DoesNotExist:
                        annotation = 'Genome not found'
                    except Exception:
                        annotation = '-'

                    if database_choice == 'bcvhp_crispr':
                        sequence_name = values[1].split('|')[1] if '|' in values[1] and len(values[1].split('|')) > 1 else '-'
                        sequence_detail = values[1].split('|')[2] if '|' in values[1] and len(values[1].split('|')) > 2 else '-'
                        sequence_type = 'CRISPR Spacer'
                    elif database_choice == 'bcvhp_trna':
                        sequence_name = values[1].split('|')[1] if '|' in values[1] and len(values[1].split('|')) > 1 else '-'
                        sequence_detail = values[1].split('|')[2] if '|' in values[1] and len(values[1].split('|')) > 2 else '-'
                        sequence_type = 'tRNA'
                    else:
                        sequence_name = values[1].split('|')[1] if '|' in values[1] and len(values[1].split('|')) > 1 else '-'
                        sequence_detail = values[1].split('|')[2] if '|' in values[1] and len(values[1].split('|')) > 2 else '-'
                        sequence_type = 'Unknown'

                    results.append({
                        'query_id': values[0],
                        'host_genome': host_genome_id,
                        'sequence_name': sequence_name,
                        'sequence_detail': sequence_detail,
                        'sequence_type': sequence_type,
                        'percent_identity': values[2],
                        'evalue': values[10],
                        'bit_score': values[11],
                        'annotation': annotation,
                        'gtdb_detail': gtdb_detail,
                        'identity_class': identity_class,
                        'evalue_class': evalue_class,
                        'bit_score_class': bit_score_class
                    })

            request.session['bcvhp_raw_result'] = ''.join(blast_results)
            request.session['bcvhp_parameters'] = {
                'program': blast_program,
                'evalue': evalue,
                'max_target_seqs': max_target_seqs,
                'culling_limit': culling_limit,
                'database': selected_db['name'],
                'database_description': selected_db['description'],
                'result_count': len(results),
                'original_count': len(results)
            }

            os.remove(query_filepath)
            os.remove(output_filepath)

            return render(request, 'tools/bcvhp_result.html', {
                'results': results,
                'parameters': request.session['bcvhp_parameters'],
                'result_count': len(results),
                'has_filtered': False,  
                'filtered_count': 0,
                'available_databases': available_databases
            })

        except subprocess.TimeoutExpired:
            if 'query_filepath' in locals() and os.path.exists(query_filepath):
                os.remove(query_filepath)
            if 'output_filepath' in locals() and os.path.exists(output_filepath):
                os.remove(output_filepath)
            request.session['bcvhp_timeout_sequence'] = sequence_input
            request.session['bcvhp_timeout_parameters'] = {
                'blast_program': blast_program,
                'evalue': evalue,
                'max_target_seqs': max_target_seqs,
                'culling_limit': culling_limit,
                'database': database_choice,
                'database_name': selected_db['name']
            }
            return redirect('bcvhp_timeout')

        except Exception as e:
            if 'query_filepath' in locals() and os.path.exists(query_filepath):
                os.remove(query_filepath)
            if 'output_filepath' in locals() and os.path.exists(output_filepath):
                os.remove(output_filepath)
            return render(request, 'tools/bcvhp.html', {
                'error': f'An unexpected error occurred: {str(e)}',
                'sequence': sequence_input,
                'blast_program': blast_program,
                'evalue': evalue,
                'max_target_seqs': max_target_seqs,
                'culling_limit': culling_limit,
                'database': database_choice,
                'available_databases': available_databases
            })

    return render(request, 'tools/bcvhp.html', {
        'available_databases': available_databases
    })

def bcvhp_timeout(request):

    timeout_sequence = request.session.get('bcvhp_timeout_sequence', '')
    timeout_parameters = request.session.get('bcvhp_timeout_parameters', {})
    return render(request, 'tools/bcvhp_timeout.html', {
        'timeout_sequence': timeout_sequence,
        'timeout_parameters': timeout_parameters
    })

def export_bcvhp_result(request):
    raw_result = request.session.get('bcvhp_raw_result', '')
    response = HttpResponse(raw_result, content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="bcvhp_blast_result.txt"'
    return response


