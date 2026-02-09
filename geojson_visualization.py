#!/usr/bin/env python

# import modules
import os
import json
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt

class MutateaGeoViz:
    def __init__(self, metadata_dir: str, output_dir: str = "geojson_output"):
        self.metadata_dir = metadata_dir
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        # create default coordinates database for Texas cities
        self.city_coordinates = {
            "Houston": {"lat": 29.7604, "lon": -95.3698, "region": "Region 6"},
            "Austin": {"lat": 30.2672, "lon": -97.7431, "region": "Region 7"},
            "Dallas": {"lat": 32.7767, "lon": -96.7970, "region": "Region 3"},
        }

    # load the processed wastewater metadata for given pathogen
    def load_metadata(self, pathogen: str) -> pd.DataFrame:
        metadata_file = os.path.join(self.metadata_dir, f"metadata_wastewater_combined.csv")
        if not os.path.exists(metadata_file):
            raise FileNotFoundError(f"Metadata file not found: {metadata_file}")
        return pd.read_csv(metadata_file)

    # create GeoJSON visualizations for given pathogen
    def create_geojson_visualizations(self, pathogen: str) -> None:
        df = self.load_metadata(pathogen)

        # add coordinates to the dataframe
        df_with_coords = self.add_coordinates(df)
        
        # create GeoJSON file
        geojson_file = f"{pathogen}_samples.geojson"
        self.create_geojson(df_with_coords, geojson_file)
        
        # create interactive HTML map
        self.create_html_map(geojson_file, pathogen)
        
        print(f"Creating visualizations for {pathogen}")
    
    def add_coordinates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add latitude and longitude coordinates based on city names"""
        df = df.copy()
        df['latitude'] = None
        df['longitude'] = None
        
        for idx, row in df.iterrows():
            city = row['City']
            if city in self.city_coordinates:
                df.at[idx, 'latitude'] = self.city_coordinates[city]['lat']
                df.at[idx, 'longitude'] = self.city_coordinates[city]['lon']
        
        return df.dropna(subset=['latitude', 'longitude'])
    
    # load coordinates for US states, counties, and world countries
    def load_geographic_data(self):    
        # counties
        try:
            self.us_counties = gpd.read_file('https://www2.census.gov/geo/tiger/GENZ2021/shp/cb_2021_us_county_20m.zip')
            self.texas_counties = self.us_counties[self.us_counties['STATEFP'] == '48']
            print("Loaded US counties data")
        except Exception as e:
            print(f"Could not load US counties data: {e}")
            self.us_counties = None
            self.texas_counties = None

        # states
        try:
            self.us_states = gpd.read_file('https://www2.census.gov/geo/tiger/GENZ2021/shp/cb_2021_us_state_20m.zip')
            print("Loaded US states data")
        except Exception as e:
            print(f"Could not load US states data: {e}")
            self.us_states = None
        
        # countries
        try:
            self.world_countries = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
            print("Loaded world countries data")
        except Exception as e:
            print(f"Could not load world countries data: {e}")
            self.world_countries = None

    def create_geojson(self, df: pd.DataFrame, output_file: str) -> None:
        """Create GeoJSON file from metadata DataFrame"""
        features = []
        
        for idx, row in df.iterrows():
            geometry = Point(row['longitude'], row['latitude'])
            properties = {
                'sample_id': row['Sample_ID'],
                'city': row['City'],
                'region': row['Region'],
                'date': str(row['Date']) if pd.notna(row['Date']) else None
            }
            
            feature = {
                'type': 'Feature',
                'geometry': geometry.__geo_interface__,
                'properties': properties
            }
            features.append(feature)
        
        geojson = {'type': 'FeatureCollection', 'features': features}
        
        output_path = os.path.join(self.output_dir, output_file)
        with open(output_path, 'w') as f:
            json.dump(geojson, f, indent=2)
        
        print(f"GeoJSON saved to: {output_path}")

    def load_geojson_data(self, geojson_file: str) -> dict:
        """Load GeoJSON data from file"""
        with open(geojson_file, 'r') as f:
            return json.load(f)
    
    def create_html_map(self, geojson_file: str, pathogen: str) -> None:
        """Create interactive HTML map from GeoJSON using Leaflet"""
        # Get the full path to the GeoJSON file
        geojson_path = os.path.join(self.output_dir, geojson_file)
        
        # Load the GeoJSON data
        geojson_data = json.dumps(self.load_geojson_data(geojson_path))
        
        html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>{pathogen} Sample Locations</title>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    
    <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
    <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
    
    <style>
        #map {{ height: 600px; width: 100%; }}
        .info {{ padding: 6px 8px; font: 14px/16px Arial, Helvetica, sans-serif; background: white; background: rgba(255,255,255,0.8); box-shadow: 0 0 15px rgba(0,0,0,0.2); border-radius: 5px; }}
        .legend {{ line-height: 18px; color: #555; }}
        .legend i {{ width: 18px; height: 18px; float: left; margin-right: 8px; opacity: 0.7; }}
    </style>
</head>
<body>
    <div id="map"></div>
    
    <script>
        // Embed GeoJSON data directly in the HTML
        var geojsonData = {geojson_data};
        
        // Initialize map centered on Texas
        var map = L.map('map').setView([31.0, -99.0], 6);
        
        // Add tile layer
        L.tileLayer('https://{{s}}.tile.openstreetmap.org/{{z}}/{{x}}/{{y}}.png', {{
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
        }}).addTo(map);
        
        // Color scheme for regions
        var colors = {{
            "Region 6": "#e41a1c",
            "Region 7": "#377eb8", 
            "Region 3": "#4daf4a",
            "Region 8": "#984ea3",
            "Region 1": "#ff7f00",
            "Region 2_3": "#ffff33",
            "Region 4_5N": "#a65628",
            "Region 9_10": "#f781bf",
            "Region 11": "#999999"
        }};
        
        // Function to get color based on region
        function getColor(region) {{
            return colors[region] || "#cccccc";
        }}
        
        // Add GeoJSON layer
        L.geoJSON(geojsonData, {{
            pointToLayer: function (feature, latlng) {{
                return L.circleMarker(latlng, {{
                    radius: 8,
                    fillColor: getColor(feature.properties.region),
                    color: "#000",
                    weight: 1,
                    opacity: 1,
                    fillOpacity: 0.8
                }});
            }},
            onEachFeature: function (feature, layer) {{
                var popupContent = "<b>Sample ID:</b> " + feature.properties.sample_id + 
                                "<br><b>City:</b> " + feature.properties.city +
                                "<br><b>Region:</b> " + feature.properties.region +
                                "<br><b>Date:</b> " + feature.properties.date;
                layer.bindPopup(popupContent);
            }}
        }}).addTo(map);
        
        // Add legend
        var legend = L.control({{position: 'bottomright'}});
        legend.onAdd = function (map) {{
            var div = L.DomUtil.create('div', 'info legend');
            div.innerHTML = '<h4>{pathogen} Regions</h4>';
            for (var region in colors) {{
                div.innerHTML +=
                    '<i style="background:' + getColor(region) + '"></i> ' + region + '<br>';
            }}
            return div;
        }};
        legend.addTo(map);
    </script>
</body>
</html>
"""
        
        # Format the string with pathogen and geojson data
        html_content = html_content.format(pathogen=pathogen, geojson_data=geojson_data)
        
        html_file = os.path.join(self.output_dir, f"{pathogen}_map.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        print(f"Interactive HTML map saved to: {html_file}")

# GeoJSON CLI (allows it to be called as module now, standalone function)
def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Create GeoJSON visualizations for mutatea outputs")
    parser.add_argument("-m", "--metadata_dir", required=True, help="Directory containing processed metadata files")
    parser.add_argument("-p", "--pathogen", required=True, help="Pathogen name (e.g., ecoli, H1N1)")
    parser.add_argument("-o", "--output_dir", default="geojson_output", help="Output directory for visualization files")
    
    args = parser.parse_args()
    
    # Create visualizations
    viz = MutateaGeoViz(args.metadata_dir, args.output_dir)
    viz.create_geojson_visualizations(args.pathogen)

if __name__ == "__main__":
    main()