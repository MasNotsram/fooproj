function OSmappingOL3(args) {
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
    var osAttribution = [new ol.Attribution({
            html: '&copy; <a href="http://www.ordnancesurvey.co.uk/">Ordnance Survey</a>'
        })];

    (function __construct(object, args) {

    }(this, args));

    this.addOSMapLayer = function (map, layerCreatedCallback, errorCallback) {
        olmap = map;
        getOSLayerCapabilities(function () {
            osLayer = new ol.layer.Tile({
                visible: true,
                source: osMapTypes.Leisure.layerSource
            });
            olmap.addLayer(osLayer);
            olmap.getView().on("change:resolution", function () {
                zoomChanged();
            });
            layerCreatedCallback(osLayer);
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
                var parser = new ol.format.WMTSCapabilities();
                var result = parser.read(data);
                $.each(osMapTypes, function (i, typeInfo) {
                    typeInfo.sourceOptions = ol.source.WMTS.optionsFromCapabilities(result, {
                        layer: typeInfo.layerName,
                        matrixSet: typeInfo.matrixSet
                    });
                    typeInfo.sourceOptions.attributions = osAttribution;
                    typeInfo.layerSource = new ol.source.WMTS(typeInfo.sourceOptions);
                });
                successCallback();
            },
            "error": function (a, b, c) {
                errorCallback(a, b, c);
            }
        });
    }

    function zoomChanged() {
        var zoom = olmap.getView().getZoom();
        osLayer.setSource(zoom < 17 ? osMapTypes.Leisure.layerSource : osMapTypes.Road.layerSource);
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