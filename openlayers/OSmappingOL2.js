function OSmappingOL2(args) {
    var olmap, osLayer;
    var osAuditSent = false;

    var osMapTypes = {
        'Leisure': {
            layerName: 'Leisure 27700',
            matrixSet: 'EPSG:27700'
        },
        'Road': {
            layerName: 'Road 27700',
            matrixSet: 'EPSG:27700'
        }
    };
    var osAttribution = "<div class='olOSAttribution'>"
            + "<div><img style='width:75px;' src='/static/openlayers/os_logo.png'></div>"
            + " &copy; Crown Copyright &amp; Database Right 2012 Ordnance Survey. All rights reserved."
            + "</div>";

    (function __construct(object, args) {

    }(this, args));

    this.addOSMapLayer = function (map, layerCreatedCallback, errorCallback) {
        olmap = map;
        getOSLayerCapabilities(function (newOSLayer) {
            osLayer = newOSLayer;
            map.addLayer(osLayer);

            layerCreatedCallback(newOSLayer);
        }, function (a, b, c) {
            errorCallback(a, b, c);
        });
    };

    function getOSLayerCapabilities(successCallback, errorCallback) {
        $.ajax({
            "type": "GET",
            "url": "https://api2.ordnancesurvey.co.uk/mapping_api/v1/service/wmts?Request=GetCapabilities&key=sxTAXVHzEliIeObIK7XeJNR9QR2wP508",
            "dataType": "text",
            "success": function (data) {
                var parser = new OpenLayers.Format.WMTSCapabilities();
                var result = parser.read(data);
                var layer = parser.createLayer(result, {
                    layer: osMapTypes.Leisure.layerName,
                    matrixSet: osMapTypes.Leisure.matrixSet,
                    format: "image/png",
                    maxExtent: new OpenLayers.Bounds(0, 0, 800000, 1300000), //needed to compute OS tiles correctly
                    attribution: osAttribution,
                    displayInLayerSwitcher: false,
                    moveTo: function () {
                        //Used to check which OS layer type to show, depending on zoom level
                        var osType = this.map.getZoom() < 10 ? osMapTypes.Leisure : osMapTypes.Road;
                        this.layer = osType.layerName;
                        this.matrixSet = osType.matrixSet;

                        OpenLayers.Layer.WMTS.prototype.moveTo.apply(this, arguments);
                    }
                });
                successCallback(layer);
            },
            "error": function (a, b, c) {
                errorCallback(a, b, c);
            }
        });
    }

    this.auditOSMapUsage = function (userId, survey) {
        if (!osAuditSent) {
            $.ajax({
                url: "https://app.bto.org/bto-ws/res/audit/osproaudit",
                dataType: "jsonp",
                crossDomain: true,
                data: {
                    proj_id: survey,
                    userid: userId
                },
                success: function () {
                    osAuditSent = true;
                },
                error: function () {

                }
            });
        }
    };
}