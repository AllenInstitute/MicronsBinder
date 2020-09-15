import neuroglancer as ngl

# Main URL is temporarily broken
# ngl.set_static_content_source(url="https://neuromancer-seung-import.appspot.com")
ngl.set_static_content_source(url="https://cjauthupdate2-dot-neuromancer-seung-import.appspot.com")
# ngl.set_static_content_source(url="https://graphene-v0-dot-neuromancer-seung-import.appspot.com/")

#SECRETS
IMG_CVPATH = "gs://microns_public_datasets/pinky100_v0/son_of_alignment_v15_rechunked"
SEG_GS_CVPATH = "gs://microns_public_datasets/pinky100_v185/seg"
CLF_CVPATH = "gs://neuroglancer/pinky100_v0/clefts/mip1_d2_1175k"
SEG_CVPATH = SEG_CL_CVPATH
MITO_CVPATH = "gs://neuroglancer/pinky100_v0/mito_seg_191220"
NUC_CVPATH = "gs://neuroglancer/pinky100_v0/nucleus/seg"

IMG_SOURCE = f"precomputed://{IMG_CVPATH}"
SEG_SOURCE = f"precomputed://{SEG_GS_CVPATH}"
MITO_SOURCE = f"precomputed://{MITO_CVPATH}"
NUC_SOURCE = f"precomputed://{NUC_CVPATH}"
CLF_CVPATH = f"precomputed://{CLF_CVPATH}"

#ngl.set_server_bind_address(bind_port=8889)
VIEWER = ngl.Viewer()


def init_std_view(viewer=VIEWER):
    with viewer.txn() as state:
        state.layers["img"] = ngl.ImageLayer(source=IMG_SOURCE)
        state.layers["seg"] = ngl.SegmentationLayer(source=SEG_SOURCE)
        state.layers["mito"] = ngl.SegmentationLayer(source=MITO_SOURCE)
        state.position.voxelSize = [4, 4, 40]


def add_layer(name, source=None, layertype="seg", linkedseg=None, color=None):
    with VIEWER.txn() as state:
        if layertype == "img":
            state.layers[name] = ngl.ImageLayer(source=source)
        elif layertype == "seg":
            state.layers[name] = ngl.SegmentationLayer(source=source)
        elif layertype == "ann":
            aux_kwargs = {}
            if linkedseg is not None:
                aux_kwargs["linkedSegmentationLayer"] = linkedseg
            if color is not None:
                aux_kwargs["annotationColor"] = color

            state.layers[name] = ngl.AnnotationLayer(**aux_kwargs)
        else:
            raise Exception(f"unknown layer type: {layertype}")


def rm_layer(name):
    with VIEWER.txn() as state:
        del state.layers[name]


def select_morph_segments(ids):
    """Select morphological segments"""
    if isinstance(ids, int):
        ids = [ids]

    with VIEWER.txn() as state:
        state.layers["seg"].segments = ids


def select_mitos(ids):
    """Select mitochondria"""
    if isinstance(ids, int):
        ids = [ids]

    with VIEWER.txn() as state:
        state.layers["mito"].segments = ids


def move_to_pos(pos, voxel_res=[4, 4, 40]):
    with VIEWER.txn() as state:
        state.position = ngl.SpatialPosition(
                             voxel_size=voxel_res,
                             voxelCoordinates=pos)


def extract_pts(layername="ann"):
    with VIEWER.txn() as state:
        annotations = state.layers[layername].annotations

    pts = list(tuple(ann.point) for ann in annotations)
    return pts


def display_pts(pts, name="ann", downsmp=1, linkedseg="seg", color=None):
    if not layer_exists(name):
        add_layer(name, layertype="ann", linkedseg=linkedseg, color=color)

    anns = [ngl.PointAnnotation(point=tuple(pt), id=f"{name}{i}")
            for (i, pt) in enumerate(pts[::downsmp])]

    with VIEWER.txn() as state:
        state.layers[name].annotations = anns


def display_lines(ptsA, ptsB, name="ann",
                  linkedseg="seg", segsA=None, segsB=None,
                  ids=None, viewer=VIEWER):
    if not layer_exists(name):
        add_layer(name, layertype="ann", linkedseg=linkedseg)

    if ids is None:
        ids = range(len(ptsA))

    if segsA is None or segsB is None:
        anns = [ngl.LineAnnotation(pointA=ptA, pointB=ptB, id=f"{i}")
                for (i, ptA, ptB) in zip(ids, ptsA, ptsB)]

    else:
        anns = [ngl.LineAnnotation(pointA=ptA, pointB=ptB, id=f"{i}",
                                   segments=[str(segA), str(segB)])
                for (i, ptA, ptB, segA, segB) in
                zip(ids, ptsA, ptsB, segsA, segsB)]

    with viewer.txn() as state:
        state.layers[name].annotations = anns


def white_background(viewer=VIEWER):
    with viewer.txn() as state:
        state.perspectiveViewBackgroundColor = "#FFFFFF"
        state.crossSectionBackgroundColor = "#FFFFFF"


def layer_exists(name):
    with VIEWER.txn() as state:
        return name in state.layers
