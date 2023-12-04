from galsim import PhotonOp
from galsim import PhotonOpBuilder


class SlitlessSpec(PhotonOp):
	r"""A photon operator that applies the dispersion effects of slitless-spectroscopy.
	
	The photons will need to have wavelengths defined in order to work.
	
	Resolution = r0 + r1*lambda + r2*lambda^2 + r3*lambda^3 + ...
	
	Parameters:
		base_wavelength:    Wavelength (in nm) represented by the fiducial photon positions
		resolution:         Resolving power per 2 pixels. Needs to be same length as base_wavelength
	"""
	# what parameters are tunable
	_req_params = {"base_wavelength": float, "barycenter": list}
	_opt_params = {"resolution": list}
	# do we need RA & Decl, or anything else?
	# _single_params = [
	# 	{"zenith_angle": Angle, 
	# 	"HA": Angle, 
	# 	"zenith_coord": CelestialCoord}
	# ]
	
	def __init__(self, base_wavelength=1_000, *resolution):
		# This matches the code in ChromaticAtmosphere.
		self.base_wavelength = base_wavelength
	
		self.resolution = np.array(resolution)
	
	def applyTo(self, photon_array, local_wcs=None, rng=None):
		"""Apply the slitless-spectroscopy disspersion to the photos
	
		Parameters:
			photon_array:   A `PhotonArray` to apply the operator to.
			local_wcs:      A `LocalWCS` instance defining the local WCS for the current photon
							bundle in case the operator needs this information.  [default: None]
			rng:            A random number generator is not used.
		"""
		#photon array has .x, .y, .wavelength, .coord, .time, ...
		if not photon_array.hasAllocatedWavelengths():
			raise GalSimError("SlitlessSpec requires that wavelengths be set")
		# if local_wcs is None:
		#     raise TypeError("SlitlessSpec requires a local_wcs to be provided to applyTo")
	
		w = photon_array.wavelength
		cenx = local_wcs.origin.x
		ceny = local_wcs.origin.y
	
		# Apply the wavelength-dependent scaling
		if self.alpha != 0.0:
			scale = (w / self.base_wavelength) ** self.alpha
			photon_array.x = scale * (photon_array.x - cenx) + cenx
			photon_array.y = scale * (photon_array.y - ceny) + ceny
	
		# Apply DCR
		shift_magnitude = dcr.get_refraction(w, self.zenith_angle, **self.kw)
		shift_magnitude -= self.base_refraction
		shift_magnitude *= radians / self.scale_unit
		sinp, cosp = self.parallactic_angle.sincos()
	
		du = -shift_magnitude * sinp
		dv = shift_magnitude * cosp
	
		dx = local_wcs._x(du, dv)
		dy = local_wcs._y(du, dv)
		photon_array.x += dx
		photon_array.y += dy
		# might need to change dxdz/dydz for the angle of travel through the detector.
	
	def __repr__(self):
		s = "galsim.SlitlessSpec(base_wavelength=%r, scale_unit=%r, alpha=%r, " % (
			self.base_wavelength,
			self.scale_unit,
			self.alpha,
		)
		s += ")"
		return s

class PhotonDCRBuilder(PhotonOpBuilder):
	"""Build a PhotonDCR
	"""
	# This one needs special handling for obj_coord
	def buildPhotonOp(self, config, base, logger):
		req, opt, single, takes_rng = get_cls_params(PhotonDCR)
		kwargs, safe = GetAllParams(config, base, req, opt, single)
		if 'sky_pos' in base:
			kwargs['obj_coord'] = base['sky_pos']
		return PhotonDCR(**kwargs)

RegisterPhotonOpType('SlitlessSpec', PhotonDCRBuilder())
