from poster.encode import multipart_encode
from poster.streaminghttp import register_openers
import urllib2

register_openers()
datagen, headers = multipart_encode({"file": open("test/mitotestnew.fasta"),'fileType':'sequences'})

request = urllib2.Request("http://mitomaster.mitomap.org/cgi-bin/websrvc.cgi", datagen, headers)

try:
	print urllib2.urlopen(request).read()
except urllib2.HTTPError, e:
	print "HTTP error: %d" % e.code
except urllib2.URLError, e:
	print "Network error: %s" % e.reason.args[1]
