#!/usr/bin/env ruby


require 'rubygems'
require 'trollop'

opts = Trollop::options do
  version "v0.0.1a"
  opt :list, "List of files", :required => true, :type => String
  opt :verbose, "Be verbose"
end

Trollop::die :list, "must exist"    if !File.exist?(opts[:list]) if opts[:list]

